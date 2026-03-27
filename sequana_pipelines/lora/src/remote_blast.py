"""Remote BLAST via NCBI's own BLAST URL API (bioservices.NCBIBlastAPI).

Submits one job per contig in parallel, polls until all jobs finish,
then writes a TSV matching the local blastn -outfmt 6 column layout
expected by the LORA report (BLAST_KEY in enums.py).
"""
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Tuple

from .enums import BLAST_KEY

# Maximum seconds to wait for a single job before giving up
_JOB_TIMEOUT = 1800  # 30 min


def _parse_fasta(fasta_path: Path) -> List[Tuple[str, str]]:
    """Return list of (seq_id, sequence) from a FASTA file."""
    sequences = []
    header, parts = None, []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    sequences.append((header, "".join(parts)))
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if header is not None:
        sequences.append((header, "".join(parts)))
    return sequences


def _parse_blast_xml(xml_text: str, query_id: str, max_hits: int, log_file=None) -> List[List]:
    """Parse BLAST XML (v1 or v2) and return rows matching BLAST_KEY order.

    Columns: qseqid sseqid pident length mismatch gapopen qstart qend
             sstart send evalue bitscore staxids stitle
    """
    # Strip BLAST XML v2 namespace for uniform parsing
    xml_text = xml_text.strip()
    if "xmlns" in xml_text[:300]:
        xml_text = xml_text.replace(' xmlns="http://www.ncbi.nlm.nih.gov"', "")

    # Results may be wrapped in a web page — extract the XML block
    start = xml_text.find("<?xml")
    if start > 0:
        xml_text = xml_text[start:]

    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        print(
            f"[remote_blast] WARNING: could not parse XML for {query_id}: {exc}",
            file=(log_file if log_file is not None else sys.stderr),
            flush=True,
        )
        return []

    rows = []
    for hit in root.iter("Hit"):
        sseqid = hit.findtext("Hit_accession") or hit.findtext("Hit_id") or ""
        stitle = hit.findtext("Hit_def") or ""
        staxids = hit.findtext("Hit_taxid") or ""

        # Best HSP is the first one (BLAST returns them sorted by score)
        hsp = next(iter(hit.iter("Hsp")), None)
        if hsp is None:
            continue

        align_len = int(hsp.findtext("Hsp_align-len") or 0)
        identity = int(hsp.findtext("Hsp_identity") or 0)
        gap_chars = int(hsp.findtext("Hsp_gaps") or 0)

        pident = f"{identity / align_len * 100:.2f}" if align_len else "0.00"
        mismatch = max(0, align_len - identity - gap_chars)

        # Count gap openings from alignment strings when available
        qseq = hsp.findtext("Hsp_qseq") or ""
        hseq = hsp.findtext("Hsp_hseq") or ""
        gapopen = _count_gap_openings(qseq) + _count_gap_openings(hseq) if qseq and hseq else gap_chars

        rows.append(
            [
                query_id,
                sseqid,
                pident,
                align_len,
                mismatch,
                gapopen,
                hsp.findtext("Hsp_query-from") or "",
                hsp.findtext("Hsp_query-to") or "",
                hsp.findtext("Hsp_hit-from") or "",
                hsp.findtext("Hsp_hit-to") or "",
                hsp.findtext("Hsp_evalue") or "",
                hsp.findtext("Hsp_bit-score") or "",
                staxids,
                stitle,
            ]
        )
        if len(rows) >= max_hits:
            break

    return rows


def _count_gap_openings(seq: str) -> int:
    """Count the number of gap-opening events (runs of '-') in an aligned sequence."""
    count, in_gap = 0, False
    for ch in seq:
        if ch == "-":
            if not in_gap:
                count += 1
                in_gap = True
        else:
            in_gap = False
    return count


def run_remote_blast(
    fasta_path: str,
    output_tsv: str,
    email: str,
    database: str = "nt",
    evalue: str = "1e-10",
    max_target_seqs: int = 100,
    max_seq_length: int = 100_000,
    api_key: str = "",
    entrez_query: str = "",
    log_file=None,
) -> None:
    """Submit all sequences in *fasta_path* to NCBI BLAST and write *output_tsv*.

    Parameters
    ----------
    fasta_path:      Input FASTA file (typically the seqkit_head subset).
    output_tsv:      Output path for the TSV results file.
    email:           Contact email required by NCBI's usage policy.
    database:        NCBI database (default: nt). Common values: nt, nr,
                     refseq_genomic, refseq_rna, refseq_protein, swissprot.
    evalue:          E-value threshold passed to the BLAST job.
    max_target_seqs: Maximum number of hits to retrieve per sequence.
    max_seq_length:  Sequences longer than this (in bp) are truncated to their
                     first *max_seq_length* bases before submission. NCBI's URL
                     API accepts up to ~1 Mb, but keeping queries short reduces
                     queue times significantly.
    api_key:         Optional NCBI API key. Registered users are subject to
                     higher rate limits and less aggressive CPU throttling.
                     Obtain one at https://www.ncbi.nlm.nih.gov/account/
    entrez_query:    Optional Entrez query to restrict the search to a subset of
                     the database — equivalent to the "Organism" box on the NCBI
                     BLAST web form. E.g. "Bacteria[Organism]" or
                     "Streptococcus[Organism]". Leave empty to search the full
                     database.
    log_file:        Open file object for progress messages. Defaults to stderr.
    """
    out = log_file if log_file is not None else sys.stderr

    def _log(msg: str) -> None:
        print(msg, file=out, flush=True)

    try:
        from bioservices import NCBIBlastAPI
    except ImportError:
        raise ImportError(
            "bioservices >= 1.16.0 is required for remote BLAST. " "Install it with: pip install 'bioservices>=1.16.0'"
        )

    if not email:
        raise ValueError(
            "An email address is required by NCBI's usage policy. "
            "Set blast.email in config.yaml or pass --blast-email on the CLI."
        )

    blast = NCBIBlastAPI(verbose=False)
    if api_key:
        blast.api_key = api_key
    sequences = _parse_fasta(fasta_path)

    if not sequences:
        _log(f"[remote_blast] No sequences found in {fasta_path}; writing empty output.")
        Path(output_tsv).write_text("")
        return

    # Truncate sequences that exceed our query length limit
    capped = []
    for seq_id, sequence in sequences:
        if len(sequence) > max_seq_length:
            _log(
                f"[remote_blast] {seq_id} ({len(sequence):,} bp) exceeds the {max_seq_length:,} bp "
                f"query limit — truncating to first {max_seq_length:,} bp to reduce NCBI queue time."
            )
            capped.append((seq_id, sequence[:max_seq_length]))
        else:
            capped.append((seq_id, sequence))

    if api_key:
        _log("[remote_blast] Using NCBI API key (higher rate limits).")
    if entrez_query:
        _log(f"[remote_blast] Restricting search with Entrez query: {entrez_query!r}")

    _log(f"[remote_blast] Submitting {len(capped)} sequence(s) to NCBI BLAST ({database}) sequentially…")

    # Submit and collect one sequence at a time to avoid IP-level CPU throttling.
    # Submitting all jobs in parallel causes NCBI to reject or silently drop results
    # for large queries against nt when multiple jobs compete from the same IP.
    _RETRY_DELAYS = [15, 30, 60]  # seconds between fetch retries for empty results
    # Delay between consecutive submissions to stay within NCBI's usage policy.
    # With an API key the limit is 10 req/s; without it is 3 req/s. We are
    # conservative here since BLAST jobs are expensive on NCBI's side.
    _INTER_JOB_DELAY = 5 if api_key else 15  # seconds

    all_rows: List[List] = []
    for i, (seq_id, sequence) in enumerate(capped):
        if i > 0:
            _log(f"[remote_blast]   waiting {_INTER_JOB_DELAY}s before next submission…")
            time.sleep(_INTER_JOB_DELAY)

        run_kwargs: Dict = dict(
            program="blastn",
            database=database,
            sequence=sequence,
            email=email,
            evalue=evalue,
            hitlist_size=max_target_seqs,
        )
        if entrez_query:
            # NCBI BLAST URL API parameter names are uppercase; kwargs are
            # forwarded verbatim by bioservices, so we must use the exact name.
            run_kwargs["ENTREZ_QUERY"] = entrez_query

        rid, rtoe = blast.run(**run_kwargs)
        _log(f"[remote_blast]   {seq_id} → RID {rid} (est. {rtoe}s)")

        status = blast.wait(rid, rtoe)
        if status == "READY":
            xml_text = blast.get_result(rid, format_type="XML")
            rows = _parse_blast_xml(xml_text, seq_id, max_hits=max_target_seqs, log_file=out)
            # Retry if NCBI returned READY but with an empty result set (can
            # happen transiently under load)
            for attempt, delay in enumerate(_RETRY_DELAYS, start=1):
                if rows:
                    break
                _log(f"[remote_blast]   {seq_id}: 0 hits on attempt {attempt}, " f"retrying in {delay}s (RID {rid})…")
                time.sleep(delay)
                xml_text = blast.get_result(rid, format_type="XML")
                rows = _parse_blast_xml(xml_text, seq_id, max_hits=max_target_seqs, log_file=out)
            _log(f"[remote_blast]   {seq_id}: {len(rows)} hit(s)")
            all_rows.extend(rows)
        elif status == "TIMEOUT":
            _log(f"[remote_blast]   {seq_id}: timed out waiting for NCBI — no results written for this contig.")
        else:
            _log(f"[remote_blast]   {seq_id}: no results (job status: {status})")

    with open(output_tsv, "w") as fh:
        for row in all_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    _log(f"[remote_blast] Done — {len(all_rows)} total hits written to {output_tsv}")
