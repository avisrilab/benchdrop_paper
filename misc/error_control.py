"""Same-length error control (Methods 2.5): does nanopore error change per-cell
isoform assignment accuracy once read length is held to the same model? Two simulated pools,
perfect and errored, drawn to the same read-length distribution over the 1,037 template transcripts.

Stages: check inputs, generate both cDNA pools with NanoSim, run the staircase chain on each
arm, compare assignment accuracy between arms.

Needs the staircase template, error profile and transcriptome index (staircase/run_depth_sweep.sh
or run_error_sweep.sh builds them), NanoSim at $NANOSIM_DIR, Bagpiper and its whitelist at
$BAGPIPER_BIN / $BAGPIPER_DIR, and the reads toolchain at $BENCHDROP_READS_BIN.
"""
from __future__ import annotations
import argparse
import json
import os
import re
import statistics
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from time import perf_counter


REV = Path(os.path.expandvars("$BENCHDROP_RESULTS"))
STA = REV / "staircase"
SILO = REV / "error_control"

CODE = Path(os.path.expandvars("$BENCHDROP_REPO/staircase"))
READS_BIN = Path(os.path.expandvars("$BENCHDROP_READS_BIN"))
PY = READS_BIN / "python"

NANOSIM = Path(os.path.expandvars("$NANOSIM_DIR"))
NANOSIM_PY = NANOSIM / ".pixi/envs/default/bin/python"
NANOSIM_SIM = NANOSIM / "src/simulator.py"
NANOSIM_LOCK = NANOSIM / "pixi.lock"

BAGPIPER = Path(os.path.expandvars("$BAGPIPER_BIN"))
WHITELIST = Path(os.path.expandvars("$BAGPIPER_DIR/test/whitelist/"
                                       "human_barcode_table.csv"))

TEMPLATE_FA = STA / "template.fa"
TEMPLATE_META = STA / "template_meta.tsv"
TXOME_MMI = STA / "txome.mmi"
PROFILE_DIR = STA / "nanosim/profile"
PROFILE_PREFIX = PROFILE_DIR / "k562"

# ---- pool generation (both pools are products of this script) ----
POOL_GEN_DIR = SILO / "pool_gen"          # main-draw + top-up scratch, logs, per-tx count dumps
P_POOL_DIR = POOL_GEN_DIR / "P"
E_POOL_DIR = POOL_GEN_DIR / "E"
P_POOL = P_POOL_DIR / "pool.fasta"        # FINAL merged (main + top-up) pool, arm P (--perfect)
E_POOL = E_POOL_DIR / "pool.fasta"        # FINAL merged (main + top-up) pool, arm E (empirical)

DEPTH = 8
ERROR_SCALE = 1.0
N_POOL = 200_000          # the main draw, per arm
THREADS = 8

MIN_READS_PER_TX = 8       # every one of the 1,037 template transcripts, each pool
TOPUP_TARGET = 10          # top up past the bar by a small margin, not to the bar exactly
TOPUP_MAX_ATTEMPTS = 6     # geometric request back-off cap for the wrap-reject retry loop

# Seeds are recorded so every invocation is fully specified. NanoSim reseeds each worker from OS
# entropy before forking, so `--seed` does not make the per-read draws bit-reproducible.
NANOSIM_SEED_P_MAIN = 101
NANOSIM_SEED_E_MAIN = 102
NANOSIM_SEED_TOPUP_BASE = 2000    # top-up call i (0-indexed, arm P's transcripts then arm E's)
                                  # uses seed BASE + i for its first attempt, +1 per retry

# arm -> the cDNA pool it draws from. Everything else is identical.
ARMS = {"P": P_POOL, "E": E_POOL}
ARM_DIR = {a: SILO / a for a in ARMS}

BARS = {
    # ---- input identity
    "template_records": 1037,        # template.fa records == template_meta rows
    "bagpiper_version": "0.1.0",

    # ---- s02: pool generation, both arms
    "main_pool_reads": 200_000,      # each arm's main draw, before top-up
    "min_reads_per_tx": 8,           # every one of 1,037 transcripts, in each final pool
    "p_pool_head_tail_zero": True,   # every --perfect read name carries head=0 and tail=0

    # ---- s03/s04: each arm ran end to end on the designed truth
    "n_cells": 10,
    "gene_cell_pairs": 2998,         # per-(cell, gene) rows with truth (evaluate_rung1.py)
    # zero fallback reads in both arms: with every transcript covered, the full-transcript
    # fallback in simulate_rung1.py never triggers. Checked two ways: the raw count is 0, and it
    # equals the truth reads of transcripts absent from the arm's pool.
    "fallback_reads_zero": True,
    "fallback_consistent": True,

    # ---- s05: the compare gates on structure only (truth and cells identical across arms,
    # zero genes excluded, the designed gene-cell-pair count); the accuracies are reported.
}


# ---------- stage machinery ----------

@dataclass
class StageResult:
    metrics: dict = field(default_factory=dict)
    passed: bool | None = None      # None for a report-only stage
    note: str = ""


class Stage:
    name: str = ""
    bar: str = ""

    def run(self) -> StageResult:
        raise NotImplementedError



_MPL_ENV = dict(os.environ, MPLCONFIGDIR=str(POOL_GEN_DIR / ".mplcache"))


def _fresh(*paths: Path):
    for p in paths:
        Path(p).unlink(missing_ok=True)


def _made(*paths: Path):
    for p in paths:
        if not Path(p).exists() or Path(p).stat().st_size == 0:
            raise RuntimeError(f"stage exited but produced no {p}")


def _run(cmd, log: Path, env=None, **kw) -> subprocess.CompletedProcess:
    """Run one tool with stderr appended to `log`."""
    with open(log, "ab") as fh:
        fh.write(("\n$ " + " ".join(str(c) for c in cmd) + "\n").encode())
        fh.flush()
        return subprocess.run([str(c) for c in cmd], stderr=fh, check=True,
                              env=env, **kw)


def _pool_length_fields(fa: Path) -> dict:
    """NanoSim read names: {tx}_{pos}_{aligned|perfect}_{i}_{F|R}_{head}_{mid}_{tail}.

    Returns counts and length summaries of the head / middle (aligned) / tail fields, so the
    two pools' length models can be compared on the ALIGNED segment, which is what the perfect
    pool emits alone."""
    n = 0
    heads, mids, tails = [], [], []
    with open(fa) as fh:
        for ln in fh:
            if ln[0] != ">":
                continue
            p = ln[1:].split("_")
            n += 1
            heads.append(int(p[5]))
            mids.append(int(p[6]))
            tails.append(int(p[7]))
    return {"reads": n,
            "median_mid": statistics.median(mids), "mean_mid": round(statistics.fmean(mids), 1),
            "mean_head": round(statistics.fmean(heads), 1),
            "mean_tail": round(statistics.fmean(tails), 1),
            "frac_head_zero": round(sum(h == 0 for h in heads) / n, 4),
            "frac_tail_zero": round(sum(t == 0 for t in tails) / n, 4),
            "median_total": statistics.median(h + m + t for h, m, t in zip(heads, mids, tails))}


def _pool_transcripts(fa: Path) -> set:
    return {ln[1:].split("_")[0] for ln in open(fa) if ln[0] == ">"}


def _fasta_counts(fa: Path) -> Counter:
    c = Counter()
    with open(fa) as fh:
        for ln in fh:
            if ln[0] == ">":
                c[ln[1:].split("_")[0]] += 1
    return c


def _iter_fasta(fa: Path):
    if not fa.exists():
        return
    hdr, buf = None, []
    with open(fa) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if ln.startswith(">"):
                if hdr is not None:
                    yield hdr, "".join(buf)
                hdr, buf = ln, []
            else:
                buf.append(ln)
    if hdr is not None:
        yield hdr, "".join(buf)


def _tx_seqs() -> dict:
    seqs, cur, buf = {}, None, []
    for ln in open(TEMPLATE_FA):
        if ln.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf)
            cur, buf = ln[1:].strip(), []
        else:
            buf.append(ln.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def _template_genes() -> dict:
    """transcript -> gene from template_meta.tsv."""
    out = {}
    with open(TEMPLATE_META) as fh:
        assert fh.readline().startswith("transcript_id")
        for ln in fh:
            t, g, _i, _l = ln.rstrip("\n").split("\t")
            out[t] = g
    return out


# ---------- s01: are the inputs all present and what they should be? ----------

class S01Inputs(Stage):
    name, bar = "s01_inputs", ("template records == meta rows == 1037; profile, whitelist, txome "
                               "index present; bagpiper 0.1.0; minimap2, samtools, NanoSim python "
                               "+ simulator.py present")
    asks = "Are the staircase inputs and the toolchain present?"

    def run(self) -> StageResult:
        SILO.mkdir(parents=True, exist_ok=True)
        n_fa = sum(1 for ln in open(TEMPLATE_FA) if ln.startswith(">"))
        n_meta = sum(1 for ln in open(TEMPLATE_META)) - 1
        ver = subprocess.run([str(BAGPIPER), "--version"], capture_output=True,
                             text=True).stdout.strip()
        tools = {t: (READS_BIN / t).exists() for t in ("minimap2", "samtools", "python")}
        present = {"profile": (PROFILE_DIR / "k562_model_profile").exists(),
                   "whitelist": WHITELIST.exists(), "txome_mmi": TXOME_MMI.exists(),
                   "nanosim_py": NANOSIM_PY.exists(), "nanosim_sim": NANOSIM_SIM.exists()}
        m = {"template_fa": n_fa, "template_meta": n_meta,
             "bagpiper_version": ver, "tools": tools, "present": present}
        ok = (n_fa == n_meta == BARS["template_records"]
              and ver.split()[-1] == BARS["bagpiper_version"]
              and all(tools.values()) and all(present.values()))
        return StageResult(metrics=m, passed=ok,
                           note=f"{n_fa} template tx, {ver}")


# ---------- s02: both NanoSim pools, main draw + per-transcript top-up ----------

def _gen_main_pool(arm: str, outdir: Path, seed: int, log: Path) -> Path:
    """The main draw for one arm (same profile, same N, same subcommand for both)."""
    outdir.mkdir(parents=True, exist_ok=True)
    out_prefix = outdir / "main"
    fa = outdir / "main_aligned_reads.fasta"
    _fresh(fa, outdir / "main_aligned_error_profile", outdir / "main_unaligned_reads.fasta")
    cmd = [NANOSIM_PY, NANOSIM_SIM, "genome", "-rg", TEMPLATE_FA, "-c", PROFILE_PREFIX,
           "-o", out_prefix, "-n", N_POOL, "--seed", seed, "-t", THREADS]
    if arm == "P":
        cmd += ["--perfect", "-s", "1"]
    _run(cmd, log, env=_MPL_ENV, cwd=outdir)
    _made(fa)
    return fa


def _topup_transcript(arm: str, t: str, seq: str, need: int, seed_start: int, outdir: Path,
                      log: Path) -> tuple[list, int, int]:
    """Top up transcript `t` to >= `need` valid (non-wrapped) reads.

    A single-transcript reference plus a tight -max can HANG (NanoSim's errored path can inflate
    the reference length consumed past the sampled draw, and the linear-genome extract_read retries
    a fixed too-long length forever against a one-contig genome), fixed with -dna_type circular,
    which always finds a placement. Circular wrapping is then common (measured ~58-76% on a 120 bp
    transcript) and a wrapped read is a seam chimera (3' spliced to 5'), so every read is checked
    against its own header (`{tx}_{pos}_..._{head}_{mid}_{tail}`) and kept only if
    pos + mid <= tx_len; the request size backs off geometrically until enough survive.
    """
    tx_len = len(seq)
    tx_fa = outdir / f"_topup_{t}.fa"
    tx_fa.write_text(f">{t}\n{seq}\n")
    kept: list = []
    wrap_rejected = 0
    n_request = max(10 * need, 150)
    seed = seed_start
    for attempt in range(TOPUP_MAX_ATTEMPTS):
        out_prefix = outdir / f"_topup_{t}_try{attempt}"
        fa = Path(str(out_prefix) + "_aligned_reads.fasta")
        cmd = [NANOSIM_PY, NANOSIM_SIM, "genome", "-rg", tx_fa, "-c", PROFILE_PREFIX,
               "-o", out_prefix, "-n", n_request, "--seed", seed, "-max", tx_len,
               "-dna_type", "circular", "-t", "1"]
        if arm == "P":
            cmd += ["--perfect", "-s", "1"]
        _run(cmd, log, env=_MPL_ENV, cwd=outdir, stdout=subprocess.DEVNULL)
        seed += 1
        for hdr, rseq in _iter_fasta(fa):
            p = hdr[1:].split("_")
            pos, mid = int(p[1]), int(p[6])
            if pos + mid <= tx_len:
                kept.append((hdr, rseq))
            else:
                wrap_rejected += 1
        for suffix in ("_aligned_reads.fasta", "_aligned_error_profile", "_unaligned_reads.fasta"):
            Path(str(out_prefix) + suffix).unlink(missing_ok=True)
        if len(kept) >= need:
            break
        n_request *= 3
    tx_fa.unlink(missing_ok=True)
    return kept, seed - seed_start, wrap_rejected


def _build_pool(arm: str, main_seed: int, topup_seed_start: int, outdir: Path) -> dict:
    log = outdir / "nanosim.log"
    _fresh(log)
    main_fa = _gen_main_pool(arm, outdir, main_seed, log)
    main_count = sum(1 for _ in open(main_fa) if _.startswith(">"))
    counts = _fasta_counts(main_fa)
    tx_seqs = _tx_seqs()
    tx_gene = _template_genes()
    deficient = [t for t in tx_gene if counts.get(t, 0) < MIN_READS_PER_TX]

    topup_metrics, topup_records = {}, []
    seed_cursor = topup_seed_start
    for t in deficient:
        need = TOPUP_TARGET - counts.get(t, 0)
        kept, used, wrap_rejected = _topup_transcript(arm, t, tx_seqs[t], need, seed_cursor,
                                                       outdir, log)
        seed_cursor += used
        topup_records.extend(kept)
        topup_metrics[t] = {"gene_id": tx_gene[t], "tx_len": len(tx_seqs[t]),
                            "main_reads": counts.get(t, 0), "topup_reads": len(kept),
                            "wrap_rejected": wrap_rejected,
                            "final_reads": counts.get(t, 0) + len(kept)}

    final_fa = ARMS[arm]
    final_fa.parent.mkdir(parents=True, exist_ok=True)
    with open(final_fa, "w") as fo:
        with open(main_fa) as fi:
            fo.write(fi.read())
        for hdr, seq in topup_records:
            fo.write(hdr + "\n" + seq + "\n")

    final_counts = _fasta_counts(final_fa)
    min_reads = min((final_counts.get(t, 0) for t in tx_gene), default=0)
    return {"main_reads": main_count, "deficient_tx": len(deficient),
            "topup_reads_total": len(topup_records),
            "wrap_rejected_total": sum(v["wrap_rejected"] for v in topup_metrics.values()),
            "seeds_used": seed_cursor - topup_seed_start,
            "min_reads_per_tx": min_reads, "final_reads": sum(final_counts.values()),
            "per_tx_topup": topup_metrics}


class S02PoolGen(Stage):
    name, bar = "s02_pool_gen", ("both arms' main NanoSim draw == 200,000 reads; every one of "
                                 "1,037 template transcripts has >= 8 reads in EACH final pool "
                                 "after per-transcript top-up; P pool head == tail == 0 throughout")
    asks = ("What do the two cDNA pools look like once every template transcript is guaranteed "
            "coverage, and how does their aligned-length model compare?")

    def run(self) -> StageResult:
        SILO.mkdir(parents=True, exist_ok=True)
        (POOL_GEN_DIR / ".mplcache").mkdir(parents=True, exist_ok=True)
        r = {}
        r["P"] = _build_pool("P", NANOSIM_SEED_P_MAIN, NANOSIM_SEED_TOPUP_BASE, P_POOL_DIR)
        # arm E's top-up seeds start well clear of arm P's range.
        r["E"] = _build_pool("E", NANOSIM_SEED_E_MAIN, NANOSIM_SEED_TOPUP_BASE + 5000,
                             E_POOL_DIR)

        pf, ef = _pool_length_fields(P_POOL), _pool_length_fields(E_POOL)
        json.dump({"P": pf, "E": ef}, open(SILO / "pool_lengths.json", "w"), indent=1)
        json.dump({a: {k: v for k, v in res.items() if k != "per_tx_topup"} for a, res in r.items()},
                  open(SILO / "pool_gen_summary.json", "w"), indent=1)
        json.dump({a: res["per_tx_topup"] for a, res in r.items()},
                  open(SILO / "pool_gen_topup_detail.json", "w"), indent=1)

        ok = (r["P"]["main_reads"] == BARS["main_pool_reads"]
              and r["E"]["main_reads"] == BARS["main_pool_reads"]
              and r["P"]["min_reads_per_tx"] >= BARS["min_reads_per_tx"]
              and r["E"]["min_reads_per_tx"] >= BARS["min_reads_per_tx"]
              and (pf["frac_head_zero"] == 1.0 and pf["frac_tail_zero"] == 1.0)
              == BARS["p_pool_head_tail_zero"])
        m = {"P": {k: v for k, v in r["P"].items() if k != "per_tx_topup"},
             "E": {k: v for k, v in r["E"].items() if k != "per_tx_topup"},
             "pool_lengths": {"P": pf, "E": ef}}
        return StageResult(metrics=m, passed=ok,
                           note=(f"P main {r['P']['main_reads']:,} + topup "
                                 f"{r['P']['topup_reads_total']:,} over {r['P']['deficient_tx']} tx "
                                 f"(min {r['P']['min_reads_per_tx']}); "
                                 f"E main {r['E']['main_reads']:,} + topup "
                                 f"{r['E']['topup_reads_total']:,} over {r['E']['deficient_tx']} tx "
                                 f"(min {r['E']['min_reads_per_tx']})"))


# ---------- s03 / s04: one arm end to end ----------

class ArmStage(Stage):
    """simulate -> barcode -> align -> count -> evaluate -> identifiability, in the arm's own dir."""
    name = "arm"
    bar = ("fallback reads == 0 (every transcript covered in both pools); 10 cells; 2,998 "
           "gene-cell pairs with truth scored; every product present")

    def __init__(self, arm: str, idx: int):
        self.arm = arm
        self.name = f"s{idx:02d}_arm_{arm}"
        self.asks = (f"Arm {arm} ({'no cDNA error' if arm == 'P' else 'empirical cDNA error'}): "
                     "what are assignment accuracy, sensitivity and fraction error at depth 8?")

    def run(self) -> StageResult:
        d = ARM_DIR[self.arm]
        d.mkdir(parents=True, exist_ok=True)
        (d / "out").mkdir(exist_ok=True)
        log = d / "arm.log"
        _fresh(log)
        for src in (TEMPLATE_FA, TEMPLATE_META, TXOME_MMI):
            dst = d / src.name
            if dst.is_symlink() or dst.exists():
                dst.unlink()
            dst.symlink_to(src)
        env = dict(os.environ, STAIRCASE_DIR=str(d), NANOSIM_POOL_FA=str(ARMS[self.arm]),
                   STAIRCASE_WL=str(WHITELIST), PATH=f"{READS_BIN}:{os.environ['PATH']}")
        sim, bc_dir, bam, counts = (d / "sim.rung1.fastq.gz", d / "out/barcode",
                                    d / "aligned.bam", d / "counts")
        _fresh(sim, d / "truth.rung1.tsv", d / "cells.rung1.tsv", bc_dir / "passed.bcd.nanopore.fa.gz",
               bam, d / "percellgene.rung1.tsv", d / "percelltx.rung1.tsv", d / "sweep_results.tsv")
        for f in counts.glob("*"):
            f.unlink()

        r = subprocess.run([str(PY), str(CODE / "simulate_rung1.py"), str(ERROR_SCALE),
                            str(DEPTH)], env=env, capture_output=True, text=True, check=True)
        open(log, "a").write(r.stdout + r.stderr)
        _made(sim, d / "truth.rung1.tsv", d / "cells.rung1.tsv")
        fallback = int(re.search(r"fallback reads[^:]*: ([\d,]+)", r.stderr).group(1).replace(",", ""))
        pool_line = re.search(r"NanoSim pool: ([\d,]+) reads over (\d+) transcripts", r.stderr)
        n_reads = int(re.search(r"reads ([\d,]+),", r.stdout).group(1).replace(",", ""))

        r = subprocess.run([str(BAGPIPER), "barcode", "--r1", str(sim), "--whitelist",
                            str(WHITELIST), "--nanopore", "--output", str(bc_dir),
                            "--threads", str(THREADS)], env=env, capture_output=True,
                           text=True, check=True)
        open(log, "a").write(r.stdout + r.stderr)
        _made(bc_dir / "passed.bcd.nanopore.fa.gz")
        bc = re.search(r"total (\d+)\s+matched (\d+)\s+small (\d+)\s+ambiguous (\d+)\s+mismatch (\d+)",
                       r.stdout + r.stderr)
        bc = dict(zip(("total", "matched", "small", "ambiguous", "mismatch"), map(int, bc.groups())))

        with open(log, "ab") as fh:
            aln = subprocess.Popen([str(READS_BIN / "minimap2"), "-ax", "map-ont", "--for-only",
                                    "-N", "200", "-p", "0.9", "-t", str(THREADS),
                                    str(TXOME_MMI), str(bc_dir / "passed.bcd.nanopore.fa.gz")],
                                   stdout=subprocess.PIPE, stderr=fh)
            subprocess.run([str(READS_BIN / "samtools"), "sort", "-n", "-@", "4", "-o", str(bam),
                            "-"], stdin=aln.stdout, stderr=fh, check=True)
            aln.wait()
            if aln.returncode:
                raise RuntimeError("minimap2 failed")
        _made(bam)
        mapped = int(subprocess.run([str(READS_BIN / "samtools"), "view", "-c", "-F", "4", str(bam)],
                                    capture_output=True, text=True, check=True).stdout.strip())

        _run([BAGPIPER, "count", "--b1", bam, "-o", counts, "--threads", THREADS], log,
             stdout=subprocess.DEVNULL)
        _made(counts / "matrix.mtx.gz")

        r = subprocess.run([str(PY), str(CODE / "evaluate_rung1.py")], env=env,
                           capture_output=True, text=True, check=True)
        open(log, "a").write(r.stdout + r.stderr)
        open(d / "evaluate.txt", "w").write(r.stdout)
        _made(d / "percellgene.rung1.tsv", d / "percelltx.rung1.tsv")
        ev = {}
        mm = re.search(r"precision ([\d.]+)\s+recall ([\d.]+)\s+\(TP (\d+) FP (\d+) FN (\d+)\)", r.stdout)
        ev.update(precision=float(mm.group(1)), recall=float(mm.group(2)),
                  TP=int(mm.group(3)), FP=int(mm.group(4)), FN=int(mm.group(5)))
        mm = re.search(r"MAE: mean ([\d.]+)\s+median ([\d.]+)", r.stdout)
        ev.update(mae_mean=float(mm.group(1)), mae_median=float(mm.group(2)))
        rows = [ln.rstrip("\n").split("\t") for ln in open(d / "percellgene.rung1.tsv")][1:]
        dom_ok = sum(int(r_[5]) for r_ in rows)
        ev.update(gene_cell_pairs=len(rows), dom_correct=dom_ok,
                  assignment_accuracy=round(dom_ok / len(rows), 4))
        n_cells = sum(1 for _ in open(d / "cells.rung1.tsv")) - 1

        r = subprocess.run([str(PY), str(CODE / "identifiability.py"), self.arm], env=env,
                           capture_output=True, text=True, check=True)
        open(d / "identifiability.txt", "w").write(r.stdout)
        open(log, "a").write(r.stderr)
        _made(d / "sweep_results.tsv", d / "cellgene_uniq.tsv")

        # Every transcript is covered in both pools (s02's bar), so `missing` should be empty
        # and fallback exactly 0.
        tx_gene = _template_genes()
        missing = sorted(set(tx_gene) - _pool_transcripts(ARMS[self.arm]))
        truth_missing = 0
        with open(d / "truth.rung1.tsv") as fh:
            fh.readline()
            for ln in fh:
                _cid, t, n = ln.rstrip("\n").split("\t")
                if t in missing:
                    truth_missing += int(n)
        with open(d / "pool_missing.tsv", "w") as fo:
            fo.write("transcript_id\tgene_id\n")
            for t in missing:
                fo.write(f"{t}\t{tx_gene[t]}\n")
        fallback_consistent = (fallback == truth_missing)
        fallback_zero = (fallback == 0)

        m = {"reads": n_reads, "pool": (pool_line.group(0) if pool_line else ""),
             "fallback_reads": fallback, "fallback_zero": fallback_zero,
             "pool_missing_tx": len(missing),
             "pool_missing_genes": len({tx_gene[t] for t in missing}),
             "truth_reads_of_missing_tx": truth_missing, "fallback_consistent": fallback_consistent,
             "barcode": bc, "mapped": mapped, "n_cells": n_cells, **ev}
        json.dump(m, open(d / "arm_metrics.json", "w"), indent=1)
        ok = (fallback_zero == BARS["fallback_reads_zero"]
              and fallback_consistent == BARS["fallback_consistent"]
              and n_cells == BARS["n_cells"] and len(rows) == BARS["gene_cell_pairs"])
        return StageResult(metrics=m, passed=ok,
                           note=(f"{n_reads:,} reads, {bc['matched']:,} barcoded, "
                                 f"accuracy {ev['assignment_accuracy']:.3f}, recall {ev['recall']:.3f}, "
                                 f"MAE {ev['mae_median']:.3f}, fallback {fallback} "
                                 f"(pool lacks {len(missing)} tx)"))


# ---------- s05: the compare ----------

def _score(d: Path, genes: set | None) -> dict:
    """Re-score one arm from its persisted per-cell tables, optionally on a gene subset.

    Same definitions as evaluate_rung1.py: detection precision/recall over (cell, template
    transcript) with est>0 / truth>0; assignment accuracy = dom_correct over gene-cell pairs with
    truth; MAE median over the same pairs."""
    TP = FP = FN = 0
    with open(d / "percelltx.rung1.tsv") as fh:
        fh.readline()
        for ln in fh:
            cid, g, t, tc, ec, _tf, _ef = ln.rstrip("\n").split("\t")
            if genes is not None and g not in genes:
                continue
            tc, ec = int(tc), float(ec)
            if tc > 0 and ec > 0:
                TP += 1
            elif tc == 0 and ec > 0:
                FP += 1
            elif tc > 0 and ec == 0:
                FN += 1
    maes, ok, n = [], 0, 0
    with open(d / "percellgene.rung1.tsv") as fh:
        fh.readline()
        for ln in fh:
            cid, g, _ni, _tr, mae, dc = ln.rstrip("\n").split("\t")
            if genes is not None and g not in genes:
                continue
            n += 1
            ok += int(dc)
            maes.append(float(mae))
    return {"gene_cell_pairs": n, "dom_correct": ok, "assignment_accuracy": round(ok / n, 4),
            "recall": round(TP / (TP + FN), 4) if TP + FN else 0.0,
            "precision": round(TP / (TP + FP), 4) if TP + FP else 0.0,
            "mae_median": round(statistics.median(maes), 4)}


def _bins(d: Path, genes: set | None) -> list:
    """Gene-level assignment accuracy by unique-sequence fraction over EXPRESSED siblings, from
    identifiability.py's persisted cellgene_uniq.tsv (designed denominator: every multi-isoform
    gene in every cell)."""
    BINS = [(0.0, 0.01), (0.01, 0.05), (0.05, 0.10), (0.10, 0.25), (0.25, 0.50), (0.50, 1.01)]
    acc = {b: [] for b in BINS}
    with open(d / "cellgene_uniq.tsv") as fh:
        fh.readline()
        for ln in fh:
            _cid, g, _ua, ue, dc = ln.rstrip("\n").split("\t")
            if genes is not None and g not in genes:
                continue
            ue = float(ue)
            for lo, hi in BINS:
                if lo <= ue < hi:
                    acc[(lo, hi)].append(int(dc))
                    break
    return [(f"{lo:.2f}-{hi:.2f}", len(v), round(sum(v) / len(v), 4) if v else float("nan"))
            for (lo, hi), v in acc.items()]


class S05Compare(Stage):
    name, bar = "s05_compare", ("truth and cells byte-identical across arms; zero genes excluded "
                                "(every template transcript covered in both pools); assignment "
                                "accuracy reported for both arms")
    asks = ("Does nanopore error at the empirical rate change per-cell isoform accuracy once read "
            "length is held to the same model?")

    def run(self) -> StageResult:
        P, E = ARM_DIR["P"], ARM_DIR["E"]
        same_truth = open(P / "truth.rung1.tsv", "rb").read() == open(E / "truth.rung1.tsv", "rb").read()
        same_cells = open(P / "cells.rung1.tsv", "rb").read() == open(E / "cells.rung1.tsv", "rb").read()
        tx_gene = _template_genes()
        excluded = set()
        for d in (P, E):
            with open(d / "pool_missing.tsv") as fh:
                fh.readline()
                excluded |= {ln.rstrip("\n").split("\t")[1] for ln in fh}
        all_genes = set(tx_gene.values())

        mp, me = _score(P, None), _score(E, None)
        armP, armE = (json.load(open(d / "arm_metrics.json")) for d in (P, E))
        delta_all = round(mp["assignment_accuracy"] - me["assignment_accuracy"], 4)

        with open(SILO / "error_control.tsv", "w") as fo:
            fo.write("gene_set\tmetric\tP_no_error\tE_empirical\n")
            for label, (a, b) in (("all 300 genes", (mp, me)),):
                fo.write(f"{label}\tgene-cell pairs scored\t{a['gene_cell_pairs']:,}\t{b['gene_cell_pairs']:,}\n")
                fo.write(f"{label}\tassignment accuracy (dominant isoform recovered)\t"
                         f"{a['assignment_accuracy']:.3f} ({a['dom_correct']:,}/{a['gene_cell_pairs']:,})\t"
                         f"{b['assignment_accuracy']:.3f} ({b['dom_correct']:,}/{b['gene_cell_pairs']:,})\n")
                fo.write(f"{label}\ttranscript recall\t{a['recall']:.3f}\t{b['recall']:.3f}\n")
                fo.write(f"{label}\tprecision\t{a['precision']:.3f}\t{b['precision']:.3f}\n")
                fo.write(f"{label}\tisoform-fraction MAE, median\t{a['mae_median']:.3f}\t{b['mae_median']:.3f}\n")
            fo.write(f"run\tbarcode matched\t{armP['barcode']['matched']:,}/{armP['barcode']['total']:,}\t"
                     f"{armE['barcode']['matched']:,}/{armE['barcode']['total']:,}\n")
            fo.write(f"run\tmapped reads\t{armP['mapped']:,}\t{armE['mapped']:,}\n")
            fo.write(f"run\tfallback reads (transcript absent from the pool)\t{armP['fallback_reads']}\t{armE['fallback_reads']}\n")
            fo.write(f"run\ttranscripts absent from the pool\t{armP['pool_missing_tx']}\t{armE['pool_missing_tx']}\n")
        bp, be = _bins(P, None), _bins(E, None)
        with open(SILO / "bins_expressed.tsv", "w") as fo:
            fo.write("bin\tn\tP_no_error\tE_empirical\n")
            for (b, n, vp), (b2, n2, ve) in zip(bp, be):
                assert b == b2 and n == n2, (b, b2, n, n2)
                fo.write(f"{b}\t{n}\t{vp:.3f}\t{ve:.3f}\n")
        _made(SILO / "error_control.tsv", SILO / "bins_expressed.tsv")
        m = {"same_truth": same_truth, "same_cells": same_cells,
             "excluded_genes": len(excluded), "all_genes": len(all_genes),
             "P": mp, "E": me, "delta_P_minus_E": delta_all,
             "bins_expressed_all": [{"bin": b, "n": n, "P": vp, "E": ve}
                                    for (b, n, vp), (_, _, ve) in zip(bp, be)]}
        # gate: structural facts only (truth/cells identity, zero exclusions, the designed
        # gene-cell-pair count); the accuracy numbers are reported.
        ok = (same_truth and same_cells and len(excluded) == 0
              and mp["gene_cell_pairs"] == BARS["gene_cell_pairs"]
              and me["gene_cell_pairs"] == BARS["gene_cell_pairs"])
        return StageResult(metrics=m, passed=ok,
                           note=(f"truth identical {same_truth}, excluded genes {len(excluded)}; "
                                 f"all 300 genes P {mp['assignment_accuracy']:.3f} vs E "
                                 f"{me['assignment_accuracy']:.3f} (delta {delta_all:+.3f})"))


STAGES = [S01Inputs(), S02PoolGen(), ArmStage("P", 3), ArmStage("E", 4), S05Compare()]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--from-stage", metavar="NAME", help="resume at this stage")
    args = ap.parse_args()
    names = [s.name for s in STAGES]
    stages = STAGES[names.index(args.from_stage):] if args.from_stage else STAGES
    SILO.mkdir(parents=True, exist_ok=True)
    metrics_path = SILO / "stage_metrics.json"
    metrics = json.loads(metrics_path.read_text()) if metrics_path.exists() else {}
    for stage in stages:
        print(f"[{stage.name}] {stage.asks}\n[{stage.name}] bar: {stage.bar}")
        t0 = perf_counter()
        result = stage.run()
        status = "ok" if result.passed is None else ("pass" if result.passed else "FAIL")
        print(f"[{stage.name}] {status} ({perf_counter() - t0:.1f}s) {result.note}".rstrip())
        metrics[stage.name] = result.metrics
        metrics_path.write_text(json.dumps(metrics, indent=2, default=str))
        if result.passed is False:
            sys.exit(f"{stage.name} failed its check; stopping")


if __name__ == "__main__":
    main()
