#!/usr/bin/env python3
"""Amortized cost model behind Table 2 (Methods 1.5): reagents, flow cell, equipment depreciation
and maintenance per run, normalized per cell. Each parameter carries a status, C (a dated public
list price or a paper-stated number), P (provisional) or F (not published here).

The 2,000-cell row needs two prices that are not published, the PIP-seq T2 kit and its library
reagents (pipseq_t2_kit, library_reagents_2k); until they are supplied, that row prints FILL.

Usage:  python3 Table2.py <out.md>
"""
import sys
from datetime import date

# status: 'C' | 'P' | 'F'
PARAMS = {
    # ONT list prices, store.nanoporetech.com/priceList.html, fetched 2026-08-13
    'flowcell_minion':   dict(v=840.00,  s='C', src='ONT store price list, 2026-08-13'),
    'flowcell_prom_pack':dict(v=4160.00, s='C', src='ONT store price list, 2026-08-13 (FLO-PRO114M pack)'),
    'prom_pack_size':    dict(v=4,       s='P', src='typical ONT pack size; confirm on the store'),
    'lsk114_kit':        dict(v=720.00,  s='C', src='ONT store price list, 2026-08-13 (SQK-LSK114, 6 rxn)'),
    'lsk114_rxns':       dict(v=6,       s='C', src='kit definition'),
    'minion_mk1d':       dict(v=3150.00, s='C', src='ONT store price list, 2026-08-13 (incl 12 mo support)'),
    'p2_solo':           dict(v=10455.00,s='P', src='ONT list value; no instrument quote held; Table 2 states it as list price'),

    # Paper-stated anchor
    'reagent_anchor_10k':dict(v=0.076,   s='C', src='capture through library prep per cell at 10k: the preprint figure with its flow cell ($1,000) removed (original cost sheet)'),
    'library_reagents_2k':dict(v=None,   s='F', src='library reagents for the 2,000-cell run; the price paid is not published'),

    # Kit price (Fluent does not publish list prices; the price paid is not published here)
    'pipseq_t2_kit':     dict(v=None,    s='F', src='PIPseq T2 3-prime v4.0 PLUS kit; the price paid is not published, ask Fluent for a quote'),
    'tenx_3prime_rxn':   dict(v=2000.00, s='P', src='commonly cited 10x 3-prime GEX per-reaction reagents; replace with a quoted list price'),
    'tenx_instrument':   dict(v=100000.00, s='P', src='approximate 10x Chromium instrument price'),

    # Amortization assumptions (stated in the table notes, tunable)
    'runs_per_year':     dict(v=24,      s='P', src='model assumption: two sequencing runs per month'),
    'amort_years':       dict(v=3,       s='P', src='model assumption: 3-year straight-line equipment depreciation'),
    'maint_train_year':  dict(v=0.00,    s='C', src='PIP-seq: no support contract; 10x row uses tenx_service_frac'),
    'tenx_service_frac': dict(v=0.30,    s='C', src='10x service contract, 30% of instrument price per year'),

    # Run-yield metrics for the per-usable-read and per-gene normalizations
    'reads_prom_run':    dict(v=145_400_000, s='C', src='PBMC PromethION, one flow cell: 145.4M single-cell long-read sequences (manuscript Results, Methods 1.3)'),
    'bc_recovery':       dict(v=0.90,    s='C', src='manuscript: >90% barcode recovery on own data (floor value used)'),
    'genes_per_cell':    dict(v=969,     s='C', src='PBMC median genes per cell over 16,221 labelled cells'),
}


def get(k):
    return PARAMS[k]['v'], PARAMS[k]['s']


def worst(*statuses):
    for s in ('F', 'P', 'C'):
        if s in statuses:
            return s
    return 'C'


def fmt(v, s, money=True):
    if v is None or s == 'F' and v is None:
        return 'FILL'
    tag = '' if s == 'C' else f' [{s}]'
    return (f'{v:,.3f}' if money and v < 1 else f'{v:,.2f}' if money else f'{v:,.0f}') + tag


def scenario(n_cells, platform):
    """Total amortized run cost and per-cell for one run at n_cells."""
    rows = {}
    # reagents: at 10k and 20k, the 0.076/cell anchor (the preprint's 0.176
    # with its flow cell removed) scaled linearly; the 2,000-cell run is the smallest PIP-seq kit
    # (T2) plus its library reagents. No kit price is printed; the row prints totals and shares.
    if n_cells <= 2000:
        kit, s_k = get('pipseq_t2_kit'); lib, s_l = get('library_reagents_2k')
        reagents, s_r = (None, 'F') if kit is None else (kit + lib, worst(s_k, s_l))
    else:
        r10k, s_r = get('reagent_anchor_10k')
        reagents = r10k * n_cells
    rows['Reagents (capture + library prep, paper anchor)'] = (reagents, s_r)

    if platform == 'MinION':
        fc, s_fc = get('flowcell_minion')
        dev, s_dev = get('minion_mk1d')
    else:
        pack, s_p = get('flowcell_prom_pack')
        size, s_n = get('prom_pack_size')
        fc, s_fc = pack / size, worst(s_p, s_n)
        dev, s_dev = get('p2_solo')
    rows['Flow cell (consumed per run)'] = (fc, s_fc)

    runs, s_runs = get('runs_per_year')
    years, s_y = get('amort_years')
    if dev is None:
        rows['Equipment (straight-line over runs)'] = (None, 'F')
        s_eq, eq = 'F', 0.0
    else:
        eq = dev / (runs * years)
        s_eq = worst(s_dev, s_runs, s_y)
        rows['Equipment (straight-line over runs)'] = (eq, s_eq)
    mt, s_mt = get('maint_train_year')
    mt_run = (mt / runs) if mt is not None else None
    rows['Maintenance + training (annual / runs)'] = (mt_run, worst(s_mt, s_runs))

    vals = [v for v, _ in rows.values()]
    total = None if any(v is None for v in vals) else sum(vals)
    s_tot = worst(*(s for _, s in rows.values()))
    return rows, total, s_tot


def scenario_tenx(n_cells):
    """10x Chromium capture plus the identical Nanopore arm (PromethION), one run at n_cells.
    Capture reagents: one 10x 3-prime reaction per 10,000 cells (kit-sized, rounded up).
    Instrument: straight-line amortization over runs x years, plus the service contract as a
    fraction of instrument price per year. The ONT arm (flow cell, LSK114 share, P2 amortization)
    is the BenchDrop-seq PromethION arm unchanged."""
    import math
    rows = {}
    rxn, s_rxn = get('tenx_3prime_rxn')
    n_rxn = max(1, math.ceil(n_cells / 10000))
    rows['10x capture reagents (reactions x list price)'] = (rxn * n_rxn, s_rxn)
    lsk, s_l = get('lsk114_kit'); lrx, s_lr = get('lsk114_rxns')
    rows['ONT library prep (LSK114 share)'] = (lsk / lrx, worst(s_l, s_lr))
    pack, s_p = get('flowcell_prom_pack'); size, s_n = get('prom_pack_size')
    rows['Flow cell (consumed per run)'] = (pack / size, worst(s_p, s_n))
    runs, s_runs = get('runs_per_year'); years, s_y = get('amort_years')
    inst, s_i = get('tenx_instrument'); p2, s_p2 = get('p2_solo')
    rows['Equipment (10x instrument + P2, straight-line over runs)'] = (
        (inst + p2) / (runs * years), worst(s_i, s_p2, s_runs, s_y))
    frac, s_f = get('tenx_service_frac')
    rows['10x service contract (fraction of instrument / runs)'] = (inst * frac / runs, worst(s_f, s_i, s_runs))
    total = sum(v for v, _ in rows.values())
    return rows, total, worst(*(s for _, s in rows.values()))


def main():
    out = [f'<!-- GENERATED by Table2.py, {date.today().isoformat()}; edit PARAMS and regenerate. -->\n']
    out.append('# Amortized cost model (Table 2)\n')
    out.append('Formula: per-cell cost at n cells = [reagents(n) + flow cell + equipment/(runs x years) + '
               '(maintenance+training)/runs] / n. Reagents at 10,000 and 20,000 cells are $0.076/cell, the '
               'preprint figure with its flow cell removed, scaled linearly; the 2,000-cell run is the '
               'smallest PIP-seq kit plus its library reagents. Status tags: [P] provisional, verify; FILL = '
               'to be supplied; untagged = confirmed (dated source in the parameter table).\n')

    out.append('## Parameters\n')
    out.append('| Parameter | Value | Status | Source |')
    out.append('|---|---|---|---|')
    for k, p in PARAMS.items():
        val = 'FILL' if p['v'] is None else f"{p['v']:,}"
        out.append(f"| {k} | {val} | {p['s']} | {p['src']} |")
    out.append('')

    out.append('## Amortized cost per run and per cell (BenchDrop-seq)\n')
    scenarios = [(2000, 'MinION'), (10000, 'PromethION'), (20000, 'PromethION')]
    out.append('| Cells | Platform | Reagents | Flow cell | Equipment | Maint+train | Total/run | Per cell |')
    out.append('|---|---|---|---|---|---|---|---|')
    for n, plat in scenarios:
        rows, total, s_tot = scenario(n, plat)
        vals = list(rows.values())
        cells = ' | '.join(fmt(v, s) for v, s in vals)
        per_cell = None if total is None else total / n
        out.append(f'| {n:,} | {plat} | {cells} | {fmt(total, s_tot)} | {fmt(per_cell, s_tot)} |')
    out.append('')
    out.append('## Shares of run cost (the relative form for the table note)\n')
    out.append('| Cells | Platform | Flow cell share | Equipment share | Reagent share |')
    out.append('|---|---|---|---|---|')
    for n, plat in scenarios:
        rows, total, _ = scenario(n, plat)
        if total is None:
            out.append(f'| {n:,} | {plat} | FILL | FILL | FILL |')
            continue
        v = {k: x for k, (x, _) in rows.items()}
        fc = v['Flow cell (consumed per run)'] / total
        eq = v['Equipment (straight-line over runs)'] / total
        rg = v['Reagents (capture + library prep, paper anchor)'] / total
        out.append(f'| {n:,} | {plat} | {fc:.0%} | {eq:.0%} | {rg:.0%} |')
    out.append('')

    out.append('## 10x Chromium plus Nanopore, same formula, PromethION arm\n')
    out.append('| Cells | Line | Cost |')
    out.append('|---|---|---|')
    for n in (10000, 20000):
        rows, total, s_tot = scenario_tenx(n)
        for k, (v, s) in rows.items():
            out.append(f'| {n:,} | {k} | {fmt(v, s)} |')
        _bd_rows, bd_total, bd_s = scenario(n, 'PromethION')
        out.append(f'| {n:,} | **Total per run** | {fmt(total, s_tot)} |')
        out.append(f'| {n:,} | **Per cell** | {fmt(total / n, s_tot)} |')
        out.append(f'| {n:,} | **Ratio to BenchDrop-seq per cell** | {total / bd_total:.2f}x [{worst(s_tot, bd_s)}] |')
    out.append('')

    out.append('## Table 2 as printed (totals and shares; no line-item dollars)\n')
    out.append('| Cells per run | Capture | Sequencer | Total per run | Per cell | Flow cell | Reagents | Equipment and service |')
    out.append('|---|---|---|---|---|---|---|---|')
    for n, plat in scenarios:
        rows, total, _ = scenario(n, plat)
        if total is None:
            out.append(f'| {n:,} | PIP-seq | {plat} | FILL | FILL | FILL | FILL | FILL |')
            continue
        v = {k: x for k, (x, _) in rows.items()}
        fc = v['Flow cell (consumed per run)'] / total; rg = v['Reagents (capture + library prep, paper anchor)'] / total
        es = (v['Equipment (straight-line over runs)'] + v['Maintenance + training (annual / runs)']) / total
        out.append(f'| {n:,} | PIP-seq | {plat} | ${total:,.0f} | ${total / n:,.2f} | {fc:.0%} | {rg:.0%} | {es:.0%} |')
    rows, total, _ = scenario_tenx(10000); v = {k: x for k, (x, _) in rows.items()}
    fc = v['Flow cell (consumed per run)'] / total
    rg = (v['10x capture reagents (reactions x list price)'] + v['ONT library prep (LSK114 share)']) / total
    es = (v['Equipment (10x instrument + P2, straight-line over runs)'] + v['10x service contract (fraction of instrument / runs)']) / total
    out.append(f'| 10,000 | 10x Chromium | PromethION | ${total:,.0f} | ${total / 10000:,.2f} | {fc:.0%} | {rg:.0%} | {es:.0%} |')
    out.append('')
    _, t10, _ = scenario(10000, 'PromethION'); reads, _ = get('reads_prom_run'); bc, _ = get('bc_recovery'); g, _ = get('genes_per_cell')
    usable = reads * bc / 10000
    out.append(f'At 10,000 cells: {usable:,.0f} usable reads per cell, ${t10 / (reads * bc) * 1000:.4f} per 1,000 usable reads, ${t10 / 10000 / g:.4f} per gene detected at a median of {g} genes per cell.\n')

    with open(sys.argv[1], 'w') as f:
        f.write('\n'.join(out))
    print(f'wrote {sys.argv[1]}', file=sys.stderr)


if __name__ == '__main__':
    main()
