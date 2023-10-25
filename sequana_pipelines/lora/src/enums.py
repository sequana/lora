BLAST_KEY = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "staxids",
    "stitle",
]
BUSCO_KEY = ["Complete BUSCOs", "Complete and duplicated BUSCOs", "Fragmented BUSCOs", "Missing BUSCOs"]
QUAST_KEY = [
    "# contigs",
    "Largest contig",
    "Total length",
    "Mapped (%)",
    "Avg. coverage depth",
]
SET_QUAST_KEY = set(QUAST_KEY)


CHECKM_KEY = [
    "Bin Id",
    "Marker lineage",
    "# genomes",
    "# markers",
    "# marker sets",
    "0",
    "1",
    "2",
    "3",
    "4",
    "5+",
    "Completeness",
    "Contamination",
    "Strain heterogeneity",
]
