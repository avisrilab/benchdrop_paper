"""Coverage panels 2D (HBA1), 3D (CORO1B), S2B (HBA2), S2C (HBZ).
Inputs: genome.fa, genes.gtf, genes.bed.gz and the BAMs listed in BAMS. Usage: python coverage_tracks.py [GENE ...]
"""
import sys
from pathlib import Path

from integrative_transcriptomics_viewer.convenience import Configuration
from integrative_transcriptomics_viewer.export import save

READS_HORIZONTAL = dict(with_reads=True, add_track_label=False, add_reads_label=False,
                        add_coverage_label=False, vertical_layout_reads=False)
READS_VERTICAL = dict(with_reads=True, add_coverage_label=False, vertical_layout_reads=True)
COVERAGE_ONLY = dict(with_reads=False, add_coverage_label=True, vertical_layout_reads=True)

PANELS = {
    # gene: (dataset, plot options)
    "HBA1": ("K562", READS_HORIZONTAL),
    "HBA2": ("K562", READS_VERTICAL),
    "HBZ": ("K562", READS_VERTICAL),
    "CORO1B": ("PBMC", COVERAGE_ONLY),
}
BAMS = {
    "K562": {"Short": "SR.K562.sort.bam", "BenchDrop-seq": "LR.K562.genome.sort.bam"},
    "PBMC": {"Short": "SR.PBMC.sort.bam", "BenchDrop-seq": "LR.PBMC.genome.sort.bam", "Bulk": "bulk.PBMC.bam"},
}

human = Configuration(genome_fasta="genome.fa", gtf_annotation="genes.gtf", bed_annotation=["genes.bed.gz"])
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

for gene in sys.argv[1:] or PANELS:
    dataset, options = PANELS[gene]
    plot = human.plot_feature(gene, bams_dict=BAMS[dataset], tighter_track=True,
                              **{"add_track_label": gene, **options})
    save(plot, str(OUT / f"coverage_{gene}.svg"))
