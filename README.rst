<<<<<<< HEAD

.. image:: https://badge.fury.io/py/sequana-lora.svg
     :target: https://pypi.python.org/pypi/sequana_lora

.. image:: http://joss.theoj.org/papers/10.21105/joss.00352/status.svg
    :target: http://joss.theoj.org/papers/10.21105/joss.00352
    :alt: JOSS (journal of open source software) DOI

.. image:: https://github.com/sequana/lora/actions/workflows/main.yml/badge.svg
   :target: https://github.com/sequana/lora/actions/workflows/main.yaml


This is the **lora** pipeline from the `Sequana <https://sequana.readthedocs.org>`_ project

:Overview: Run assembler (Canu, flye, hifiasm) on a set of long read files
:Input: A set of BAM files from Pacbio sequencers, or FastQ files for Nanopore sequencers.
:Output: HTML reports with assemblies for each sample.
:Status: prod
:Citation: Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16), 352, JOSS DOI doi:10.21105/joss.00352
:Citation(pipeline):
    .. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.18337877.svg
        :target: https://doi.org/10.5281/zenodo.18337877
    https://www.biorxiv.org/content/10.64898/2026.01.06.697901v1

.. image:: https://raw.githubusercontent.com/sequana/lora/master/sequana_pipelines/lora/dag.png
=======
LORA — Long Read Assembly pipeline
====================================

:Overview: Assemble long reads (PacBio HiFi, PacBio subreads, Nanopore) into
           high-quality genome assemblies with optional polishing, annotation,
           and quality assessment.
:Input: BAM files from PacBio sequencers, or FastQ files from Nanopore or
        PacBio HiFi sequencers.
:Output: HTML reports with per-sample assembly statistics, coverage, BLAST
         identification, BUSCO scores, and optional annotation.
:Status: Production
:Citation: Cokelaer et al, (2017), 'Sequana': a Set of Snakemake NGS pipelines,
           Journal of Open Source Software, 2(16), 352,
           `doi:10.21105/joss.00352 <https://doi.org/10.21105/joss.00352>`_


.. image:: https://raw.githubusercontent.com/sequana/lora/master/sequana_pipelines/lora/dag.png
   :alt: Pipeline DAG
>>>>>>> d4b72d7 (Refactor pipeline rules, update README and report templates)


Installation
------------

::

    pip install sequana-lora

To upgrade an existing installation::

    pip install sequana-lora --upgrade


Quick Start
-----------

**Step 1 — prepare the working directory**::

    sequana_lora \
        --input-directory /path/to/reads \
        --data-type pacbio-hifi \
        --assembler flye \
        --genome-size 3m \
        --apptainer-prefix /path/to/containers

This creates a ``lora/`` working directory containing ``config.yaml`` and a
``lora.sh`` launch script.

**Step 2 — review the configuration** (optional but recommended)::

    cd lora
    cat config.yaml   # adjust parameters as needed

**Step 3 — run the pipeline**::

    sh lora.sh

<<<<<<< HEAD
=======
Or launch directly from step 1 with ``--execute`` (skips the review step)::

    sequana_lora ... --execute

To watch live progress in the terminal, add ``--monitor``::

    sequana_lora ... --execute --monitor

>>>>>>> d4b72d7 (Refactor pipeline rules, update README and report templates)

Required options
----------------

Three options are always required:

``--assembler``
    Assembler to use. Choices: ``flye`` (recommended for HiFi),
    ``canu``, ``hifiasm``, ``unicycler``, ``necat``, ``pecat``.

``--data-type``
    Technology and quality of the input reads:

    ============== ===============================================================
    Value          Description
    ============== ===============================================================
    pacbio-hifi    PacBio HiFi / CCS reads (Q20+)
    pacbio-raw     PacBio CLR / subreads (raw)
    pacbio-corr    PacBio corrected reads
    nano-hq        Nanopore Q20+ reads (e.g. R10.4 with SUP basecalling)
    nano-raw       Nanopore standard reads
    nano-corr      Nanopore corrected reads
    ============== ===============================================================

``--genome-size``
    Estimated genome size, e.g. ``3m`` (3 Mb), ``2.5g`` (2.5 Gb).
    Required by Flye; used by Canu for coverage reporting.


Common Examples
---------------

PacBio HiFi (recommended setup)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    sequana_lora \
        --input-directory /data/hifi \
        --data-type pacbio-hifi \
        --assembler flye \
        --genome-size 3m \
        --apptainer-prefix /shared/containers \
        --do-coverage

Nanopore (bacteria, full quality pipeline)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``--mode bacteria`` to enable sequana_coverage, prokka, busco, and checkm
in one shot::

    sequana_lora \
        --input-directory /data/nanopore \
        --data-type nano-hq \
        --assembler flye \
        --genome-size 3m \
        --apptainer-prefix /shared/containers \
        --mode bacteria \
        --busco-lineage bacteria \
        --checkm-rank genus \
        --checkm-name Streptococcus

PacBio subreads (with CCS construction)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If your input is raw PacBio BAM files (subreads), LORA can build CCS/HiFi reads
first::

    sequana_lora \
        --input-directory /data/subreads \
        --data-type pacbio-raw \
        --assembler flye \
        --genome-size 3m \
        --pacbio-build-ccs \
        --pacbio-ccs-min-passes 10 \
        --pacbio-ccs-min-rq 0.99

Or, if you have multiple BAM files per sample, provide a CSV::

    sequana_lora \
        --pacbio-input-csv samples.csv \
        --data-type pacbio-raw \
        --assembler flye \
        --genome-size 3m

The CSV format is: one row per sample with columns ``sample,file1[,file2,...]``.


Optional Steps
--------------

Coverage analysis
~~~~~~~~~~~~~~~~~

Computes depth of coverage and breadth of coverage for each contig using
sequana_coverage. Highly recommended to check assembly quality::

    --do-coverage

BUSCO completeness
~~~~~~~~~~~~~~~~~~

Assess genome completeness against a lineage-specific marker gene set::

    --busco-lineage bacteria          # auto-download bacteria lineage
    --busco-lineage streptococcales   # specific clade
    --busco-print-lineages            # list all available lineages

Prokka annotation
~~~~~~~~~~~~~~~~~

Annotate contigs (bacterial genomes)::

    --do-prokka

CheckM genome quality
~~~~~~~~~~~~~~~~~~~~~

Estimate completeness and contamination for bacterial genomes::

    --checkm-rank genus --checkm-name Streptococcus

Use an invalid ``--checkm-name`` value to get a list of valid names for a
given rank, e.g. ``--checkm-rank genus --checkm-name HELP``.

Polypolish (Illumina polishing)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Polish long-read contigs with paired-end Illumina data::

    --do-polypolish \
    --polypolish-input-directory /data/illumina \
    --polypolish-input-pattern "*.fastq.gz" \
    --polypolish-input-readtag "_R[12]_"

Circularisation
~~~~~~~~~~~~~~~

Explicit circularisation with Circlator (Flye performs this automatically)::

    --do-circlator


BLAST identification
--------------------

BLAST aligns each contig against a nucleotide database to identify the
assembled sequences. The top hits appear in the HTML report.

Local BLAST
~~~~~~~~~~~

Requires a locally installed BLAST+ and a downloaded ``nt`` database
(~270 GB). Fastest option with no network dependency::

    --blastdb /path/to/blast/databases

Remote BLAST (NCBI)
~~~~~~~~~~~~~~~~~~~

No local database required — jobs are submitted to NCBI's BLAST servers.
Enable by providing an email address::

    --blast-email your@email.com

Jobs are submitted sequentially (one contig at a time) to avoid IP-level CPU
throttling by NCBI. The default database is ``nt``; use ``--blast-remote-db``
to change it::

    --blast-email your@email.com --blast-remote-db refseq_genomic

**Restricting the search to an organism group (recommended)**

Searching all of ``nt`` for a large contig is slow and prone to NCBI throttling.
Restrict the search to a taxonomic group by editing ``config.yaml`` after
running ``sequana_lora``::

    blast:
        entrez_query: 'Bacteria[Organism]'     # all bacteria
        # entrez_query: 'Streptococcus[Organism]'   # single genus

This is equivalent to filling the "Organism" box on the NCBI BLAST web form.
It dramatically reduces search time and avoids throttling.

**NCBI API key (optional but recommended)**

Register for a free NCBI API key at https://www.ncbi.nlm.nih.gov/account/
(sign in → API Key Management). It raises the rate limit from 3 to 10
requests/second and reduces CPU throttling for large queries. Add it to
``config.yaml``::

    blast:
        api_key: 'YOUR_KEY_HERE'


HPC / SLURM cluster
-------------------

On a cluster with SLURM, pass ``--profile slurm``::

    sequana_lora \
        --input-directory /data/hifi \
        --data-type pacbio-hifi \
        --assembler flye \
        --genome-size 3m \
        --profile slurm \
        --slurm-queue fast \
        --jobs 40 \
        --apptainer-prefix /shared/containers

Per-rule memory and thread settings are controlled via the ``resources`` blocks
in ``config.yaml``.


Apptainer / Singularity (no system installs needed)
----------------------------------------------------

Every tool runs inside a pre-built container. Point ``--apptainer-prefix`` to a
shared directory so images are downloaded once and reused across projects::

    --apptainer-prefix /shared/containers

Images are downloaded automatically on first run from Zenodo. Pass extra bind
mounts with ``--apptainer-args`` if your data lives outside ``$HOME``::

    --apptainer-args "-B /data:/data"


Configuration file
------------------

After running ``sequana_lora``, a ``config.yaml`` is created in the working
directory. All pipeline parameters can be tuned there. Key sections:

- ``assembler`` — which assembler to use
- ``flye`` / ``canu`` / ``hifiasm`` — assembler-specific options
- ``fastp`` — read filtering (minimum length, etc.)
- ``blast`` — BLAST settings including ``entrez_query`` and ``api_key``
- ``busco`` / ``prokka`` / ``checkm`` — optional QC tools
- ``sequana_coverage`` — coverage analysis parameters
- ``multiqc`` — aggregated report settings

Full reference:
`config.yaml <https://raw.githubusercontent.com/sequana/sequana_lora/master/sequana_pipelines/lora/config.yaml>`_


Pipeline overview
-----------------

1. **Read filtering** — fastp removes reads below the minimum length threshold.
2. **[Optional] CCS** — build HiFi reads from PacBio subreads (ccs tool).
3. **Assembly** — Flye / Canu / Hifiasm / Unicycler / NECAT / PECAT.
4. **[Optional] Circularisation** — Circlator (or built into Flye).
5. **[Optional] Polishing** — Polypolish with paired-end Illumina reads.
6. **Contig sorting** — SeqKit sorts contigs by length (largest first).
7. **Read mapping** — Minimap2 maps reads back to contigs; Mosdepth computes coverage.
8. **[Optional] Coverage analysis** — sequana_coverage per contig.
9. **Quality assessment** — QUAST assembly statistics.
10. **[Optional] BLAST** — top hits per contig (local or remote NCBI).
11. **[Optional] BUSCO** — genome completeness.
12. **[Optional] Prokka** — genome annotation.
13. **[Optional] CheckM** — contamination and completeness for bacteria.
14. **Reports** — per-sample HTML report and a multi-sample summary.


Changelog
---------

========= ====================================================================
Version   Description
========= ====================================================================
1.1.0
1.0.0     * uniformised extension with other pipelines. fix regression on
            schema file
          * update sequana container to v0.16.5
          * add unicycler apptainer
          * add checkm module to help users choosing correct marker and name
          * replaces --pacbio and --nanopore with --data-type. pacbio is now
            decomposed into 3 sub-categories: pacbio-raw, pacbio-hifi and
            pacbio-corr
          * add bandage if assembly graph is available
          * fixed hifiasm container to use newest version
          * improved report html
          * make genome-size compulsory
          * add fastp as preprocessing tool
          * remove presets in favor of click options
          * CCS defaults to hifi. pacbio presets in config set to pacbio-hifi
          * blast removed from default; users must set blast DB themselves
          * busco lineage downloaded from the web
          * CANU preset changes: pacbio → pacbio-hifi
          * CANU-correction preset changes: pacbio → pacbio-hifi
          * FLYE preset changes: pacbio-raw → pacbio-hifi
          * remote BLAST via NCBI URL API (no local database needed)
          * entrez_query support to restrict BLAST to a taxonomic group
          * NCBI API key support for higher rate limits
0.3.0     * Use click instead of argparse
          * added multiqc / checkm / unicycler
0.2.0     * add apptainers in most rules
          * remove utils.smk to move rulegraph inside main pipeline
          * rename lora.smk into lora.rules for consistency with other
            pipelines
          * add checkm in the pipeline and HTML report
0.1.0     **First release.**
========= ====================================================================
