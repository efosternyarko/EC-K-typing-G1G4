#!/usr/bin/env python3
"""
make_v05_db.py — Build v0.5 database by replacing KL306 and KL307 with
better representative sequences identified by NCBI BLAST.

Replacements validated by running Kaptive normalised scoring against the v0.4 database:
  KL306  ESC_CC8109AA_AS (BSI) → CP099041   (Norm AS 2.0000, 31/31 CDS)
  KL307  ESC_NB5901AA_AS (BSI) → CP070103   (Norm AS 2.0000, 34/34 CDS)

KL300 replacement (CP135488) was evaluated but reverted:
  - With 1.2× extraction cap (27/30 CDS): new reference is too similar to other KX01
    variable regions — pulls many assemblies incorrectly to KL300 in normalised scoring
  - With 1.5× cap (30/30 CDS, incl. peripheral rfaG/gmd_2/KL300_11): standard Kaptive
    regresses by 54 typeable assemblies due to reference size-bias and locus fragmentation
  The old BSI KL300 representative is retained for stability.

KL303 is left unchanged — no NCBI candidate typed as KL303; the variable region
is genuinely indistinguishable from KL302/KL352 across all publicly available genomes.

Pipeline per replacement:
  1. blastn — locate reference locus in candidate chromosome → extract region
  2. blastn liftover — transfer CDS annotations from v0.4 reference by per-CDS BLAST
  3. strip conserved CDS (galF, galF_2, gnd, ugd, wza, wzb, wzc) — same as v0.4
  4. splice into v0.4 G1/G4 GenBank → write v0.5

Usage
-----
    python scripts/make_v05_db.py

Outputs (written to DB/)
------------------------
    EC-K-typing_group1and4_v0.5.gbk
    EC-K-typing_all_groups_v0.5.gbk
"""

import os
import re
import subprocess
import sys
import tempfile
from collections import OrderedDict
from datetime import date
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_DIR     = Path(__file__).resolve().parent.parent
DB_DIR       = REPO_DIR / "DB"
GENOMES_DIR  = DB_DIR / "blast_ncbi_results" / "candidate_genomes"
INPUT_G1G4   = DB_DIR / "EC-K-typing_group1and4_v0.4.gbk"
INPUT_G2G3   = DB_DIR / "EC-K-typing_group2and3_v3.0.0.gbk"
OUTPUT_G1G4  = DB_DIR / "EC-K-typing_group1and4_v0.5.gbk"
OUTPUT_ALL   = DB_DIR / "EC-K-typing_all_groups_v0.5.gbk"

# ---------------------------------------------------------------------------
# Replacements: locus → best NCBI candidate
# ---------------------------------------------------------------------------
REPLACEMENTS = {
    # KL300: CP135488 replacement tested but reverted — see notes in README.
    # With 1.2× cap (27 CDS): new reference pulls many other loci to KL300 (regression).
    # With 1.5× cap (30 CDS, incl. peripheral rfaG/gmd_2): standard Kaptive regresses
    #   (-54 typeable assemblies) due to size-bias and locus fragmentation.
    # The old BSI KL300 representative is retained for stability.
    "KL306": "CP099041",   # CP147038 wraps chromosome ORI; CP099041 is continuous
    "KL307": "CP070103",
}

# Conserved genes to strip (same set as v0.4)
CONSERVED_GENE_NAMES = {"galF", "galF_2", "gnd", "ugd", "wza", "wzb", "wzc"}

# Padding either side of BLAST match when extracting locus region
FLANK_PAD = 500


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_locus_name(rec) -> str:
    for f in rec.features:
        if f.type == "source":
            for note in f.qualifiers.get("note", []):
                m = re.search(r"K locus:\s*(\S+)", note)
                if m:
                    return m.group(1)
    return rec.name.split("_")[0]


def blast_locate_locus(ref_seq: str, chrom_fasta: Path) -> tuple[str, int, int, int]:
    """
    BLAST ref_seq against chrom_fasta.
    Returns (contig_id, start_0based, end_0based, strand +1/-1).

    Uses ALL HSPs on the best-scoring contig/strand to define locus span.
    This handles cases where the reference has internal variation that splits
    a single continuous locus into multiple BLAST HSPs.
    """
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as qf:
        qf.write(f">query\n{ref_seq}\n")
        query_fa = qf.name

    result = subprocess.run(
        ["blastn", "-query", query_fa,
         "-subject", str(chrom_fasta),
         "-outfmt", "6 sseqid sstart send sstrand bitscore pident",
         "-max_target_seqs", "1",
         "-perc_identity", "85"],
        capture_output=True, text=True,
    )
    os.unlink(query_fa)

    # Collect all HSPs; group by (contig, strand)
    from collections import defaultdict
    groups = defaultdict(list)
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts    = line.split("\t")
        contig   = parts[0]
        sstart   = int(parts[1])
        send     = int(parts[2])
        sstrand  = parts[3]
        bitscore = float(parts[4])
        groups[(contig, sstrand)].append((sstart, send, bitscore))

    if not groups:
        return None, 0, 0, 1

    # Pick the (contig, strand) with greatest total bitscore
    best_key  = max(groups, key=lambda k: sum(h[2] for h in groups[k]))
    best_hsps = groups[best_key]
    contig, sstrand = best_key

    # Filter to significant HSPs: keep only those with bitscore >= 5% of the
    # maximum bitscore. This excludes the many small spurious hits from conserved
    # short sequences scattered across the chromosome.
    max_bs    = max(h[2] for h in best_hsps)
    threshold = max(max_bs * 0.05, 1000.0)   # at least 1000, or 5% of max
    sig_hsps  = [h for h in best_hsps if h[2] >= threshold]

    # Additionally restrict to HSPs within 2× ref_len of the best HSP,
    # to exclude coincidental high-scoring distant matches.
    ref_len    = len(ref_seq)
    best_hsp   = max(sig_hsps, key=lambda h: h[2])
    best_center = (best_hsp[0] + best_hsp[1]) / 2
    max_dist    = ref_len * 2
    sig_hsps    = [h for h in sig_hsps
                   if abs((h[0] + h[1]) / 2 - best_center) <= max_dist]

    all_coords = [c for h in sig_hsps for c in (h[0], h[1])]
    raw_start  = min(all_coords) - 1   # 0-based
    raw_end    = max(all_coords)

    # Hard cap: extracted span must not exceed 1.2× ref length.
    # If the spanning HSPs cover more than that, trim inward from each end
    # keeping the primary (highest-bitscore) HSP centred.
    # Note: KL300 in CP135488 spans 63,880 bp (1.41× ref), so 1.5× would capture
    # all 30 CDS — but the 3 distal CDS (rfaG, gmd_2, KL300_11) are chromosomal
    # flanking genes outside the core locus, causing standard Kaptive locus-finding
    # to fragment the locus and flag assemblies as Untypeable. Keeping 1.2× captures
    # the 27 core variable CDS while avoiding this size-bias regression.
    max_span = int(ref_len * 1.2)
    if raw_end - raw_start > max_span:
        best_center = int((best_hsp[0] + best_hsp[1]) / 2)
        half        = max_span // 2
        raw_start   = max(0, best_center - half) - 1
        raw_end     = best_center + half

    strand  = 1 if sstrand == "plus" else -1
    return contig, raw_start, raw_end, strand


def extract_locus_from_chrom(chrom_fasta: Path, contig_id: str,
                              start_0: int, end_0: int,
                              strand: int, pad: int = FLANK_PAD) -> SeqRecord:
    """Return a SeqRecord for the extracted locus region (padded)."""
    chroms = {r.id: r for r in SeqIO.parse(str(chrom_fasta), "fasta")}
    rec    = chroms[contig_id]
    s = max(0, start_0 - pad)
    e = min(len(rec.seq), end_0 + pad)
    subseq = rec.seq[s:e]
    if strand == -1:
        subseq = subseq.reverse_complement()
    return SeqRecord(subseq, id=chrom_fasta.stem, description="")


def liftover_annotation(ref_rec: SeqRecord, new_seq: Seq) -> list:
    """
    Transfer CDS annotations from ref_rec onto new_seq by nucleotide BLAST.

    For each CDS in the reference, BLAST its nucleotide sequence against new_seq
    to find its position. Returns list of (start, end, strand, gene_name, product)
    in new_seq coordinates.

    At >98% identity the positions are essentially 1:1; any CDS that doesn't map
    is skipped (shouldn't happen for NCBI candidates at this identity level).
    """
    # Write new sequence as the BLAST subject
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as sf:
        sf.write(f">new_seq\n{str(new_seq)}\n")
        subject_fa = sf.name

    # Write all reference CDS as queries
    ref_cds = [(i, f) for i, f in enumerate(ref_rec.features) if f.type == "CDS"]
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as qf:
        for i, f in ref_cds:
            nuc = str(f.location.extract(ref_rec.seq))
            qf.write(f">cds_{i}\n{nuc}\n")
        query_fa = qf.name

    result = subprocess.run(
        ["blastn", "-query", query_fa, "-subject", subject_fa,
         "-outfmt", "6 qseqid sstart send sstrand pident",
         "-perc_identity", "85", "-max_hsps", "1", "-max_target_seqs", "1"],
        capture_output=True, text=True,
    )
    os.unlink(subject_fa)
    os.unlink(query_fa)

    # Parse: qseqid → (new_start_0, new_end_0, strand, pident)
    cds_map = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts  = line.split("\t")
        qid    = parts[0]           # cds_i
        sstart = int(parts[1])
        send   = int(parts[2])
        strd   = 1 if parts[3] == "plus" else -1
        pident = float(parts[4])
        idx    = int(qid.split("_")[1])
        if idx not in cds_map:
            start_0 = min(sstart, send) - 1
            end_0   = max(sstart, send)
            cds_map[idx] = (start_0, end_0, strd, pident)

    # Build output list in order, skipping unmapped CDS
    result_genes = []
    unmapped = 0
    for orig_i, f in ref_cds:
        gene    = f.qualifiers.get("gene",    [None])[0]
        product = f.qualifiers.get("product", ["hypothetical protein"])[0]
        if orig_i in cds_map:
            start_0, end_0, strand, _ = cds_map[orig_i]
            result_genes.append((start_0, end_0, strand, gene, product))
        else:
            unmapped += 1

    if unmapped:
        print(f"    WARNING: {unmapped} CDS could not be mapped (skipped)")

    return result_genes


def build_genbank_record(locus_name: str, sequence: Seq,
                          source_label: str, gene_list: list) -> SeqRecord:
    """Build a Kaptive-compatible GenBank record."""
    record = SeqRecord(
        Seq(str(sequence)),
        id          = f"{locus_name}_{source_label}",
        name        = f"{locus_name}_{source_label}",
        description = (
            f"Escherichia coli capsular polysaccharide synthesis "
            f"gene cluster, {locus_name}"
        ),
    )
    record.annotations.update({
        "molecule_type":      "DNA",
        "topology":           "linear",
        "data_file_division": "BCT",
        "date":               date.today().strftime("%d-%b-%Y").upper(),
        "organism":           "Escherichia coli",
        "taxonomy": [
            "Bacteria", "Pseudomonadota", "Gammaproteobacteria",
            "Enterobacterales", "Enterobacteriaceae", "Escherichia",
        ],
        "source": "Escherichia coli",
    })

    src_q = OrderedDict()
    src_q["organism"] = ["Escherichia coli"]
    src_q["mol_type"] = ["genomic DNA"]
    src_q["note"]     = [f"K locus: {locus_name}"]
    record.features.append(
        SeqFeature(FeatureLocation(0, len(sequence)), type="source", qualifiers=src_q)
    )

    gene_name_counts = {}
    for i, (start, end, strand, gene_name, product) in enumerate(gene_list):
        locus_tag = f"{locus_name}_{str(i+1).zfill(5)}"

        display_gene = gene_name
        if gene_name:
            cnt = gene_name_counts.get(gene_name, 0) + 1
            gene_name_counts[gene_name] = cnt
            if cnt > 1:
                display_gene = f"{gene_name}_{cnt}"

        q = OrderedDict()
        q["locus_tag"] = [locus_tag]
        q["gene"]      = [display_gene if display_gene else locus_tag]
        q["product"]   = [product]
        record.features.append(
            SeqFeature(FeatureLocation(start, end, strand=strand),
                       type="CDS", qualifiers=q)
        )
    return record


def strip_conserved_cds(rec: SeqRecord) -> SeqRecord:
    """Remove CDS features whose gene is in CONSERVED_GENE_NAMES."""
    rec.features = [
        f for f in rec.features
        if not (f.type == "CDS" and
                f.qualifiers.get("gene", [""])[0] in CONSERVED_GENE_NAMES)
    ]
    return rec


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 65)
    print("make_v05_db.py — Replace KL300, KL306, KL307 representatives")
    print("=" * 65)

    # Load all v0.4 G1/G4 records; collect reference records for BLASTp DB
    print(f"\n[1] Loading v0.4 G1/G4 GenBank...")
    v04_records  = list(SeqIO.parse(str(INPUT_G1G4), "genbank"))
    v04_by_name  = {get_locus_name(r): r for r in v04_records}
    print(f"    {len(v04_records)} loci loaded")

    # Process each replacement using annotation liftover (blastn of each CDS)
    new_records = {}
    for locus_name, acc in REPLACEMENTS.items():
        chrom_fasta = GENOMES_DIR / f"{locus_name}_{acc}.fasta"
        if not chrom_fasta.exists():
            print(f"\nERROR: {chrom_fasta} not found — run test_blast_candidates.py first",
                  file=sys.stderr)
            sys.exit(1)

        print(f"\n[2] Processing {locus_name} → {acc}")
        ref_rec  = v04_by_name[locus_name]
        ref_seq  = str(ref_rec.seq)
        ref_len  = len(ref_seq)
        ref_cds_count = sum(1 for f in ref_rec.features if f.type == "CDS")

        # Locate reference in chromosome
        print(f"    Locating reference ({ref_len:,} bp) in {chrom_fasta.name}...")
        contig_id, start_0, end_0, strand = blast_locate_locus(ref_seq, chrom_fasta)
        if contig_id is None:
            print(f"    ERROR: could not locate {locus_name} in {acc}", file=sys.stderr)
            sys.exit(1)
        match_len = end_0 - start_0
        print(f"    Found: {contig_id}:{start_0:,}–{end_0:,} ({match_len:,} bp, "
              f"strand {'+'if strand==1 else '-'})")

        # Extract locus region (+/- FLANK_PAD)
        locus_rec = extract_locus_from_chrom(
            chrom_fasta, contig_id, start_0, end_0, strand
        )
        extracted_len = len(locus_rec.seq)
        print(f"    Extracted: {extracted_len:,} bp (with {FLANK_PAD} bp flanking each side)")

        # Liftover annotation from v0.4 reference by blastn of each CDS
        print(f"    Lifting over {ref_cds_count} CDS annotations from v0.4 reference...")
        gene_list = liftover_annotation(ref_rec, locus_rec.seq)
        n_mapped = len(gene_list)
        n_cons   = sum(1 for g in gene_list if g[3] in CONSERVED_GENE_NAMES)
        print(f"    {n_mapped}/{ref_cds_count} CDS mapped, "
              f"{n_cons} conserved (will be stripped)")

        # Build GenBank record
        full_rec = build_genbank_record(locus_name, locus_rec.seq, acc, gene_list)

        # Strip conserved CDS (same as v0.4 make_v04_db.py)
        stripped_rec = strip_conserved_cds(full_rec)
        n_variable = sum(1 for f in stripped_rec.features if f.type == "CDS")
        print(f"    Final: {n_variable} variable CDS")

        new_records[locus_name] = stripped_rec

    # Splice new records into v0.4 G1/G4 GenBank
    print(f"\n[4] Building v0.5 G1/G4 GenBank...")
    out_g1g4 = []
    for rec in v04_records:
        lname = get_locus_name(rec)
        if lname in new_records:
            old_cds = sum(1 for f in rec.features if f.type == "CDS")
            new_cds = sum(1 for f in new_records[lname].features if f.type == "CDS")
            print(f"    {lname}: {len(rec.seq):,} bp/{old_cds} CDS  →  "
                  f"{len(new_records[lname].seq):,} bp/{new_cds} CDS  [{REPLACEMENTS[lname]}]")
            out_g1g4.append(new_records[lname])
        else:
            out_g1g4.append(rec)

    SeqIO.write(out_g1g4, str(OUTPUT_G1G4), "genbank")
    print(f"    Written: {OUTPUT_G1G4.name}  ({len(out_g1g4)} loci)")

    # Combined all-groups
    g2g3 = list(SeqIO.parse(str(INPUT_G2G3), "genbank"))
    all_recs = g2g3 + out_g1g4
    SeqIO.write(all_recs, str(OUTPUT_ALL), "genbank")
    print(f"    Written: {OUTPUT_ALL.name}  ({len(all_recs)} loci)")

    # Summary
    print("\n" + "=" * 65)
    print("SUMMARY")
    print("=" * 65)
    print(f"{'Locus':<8}  {'Old seq (bp)':>13}  {'New seq (bp)':>13}  "
          f"{'Old CDS':>8}  {'New CDS':>8}  {'Source'}")
    print("-" * 65)
    for lname, acc in REPLACEMENTS.items():
        old = v04_by_name[lname]
        new = new_records[lname]
        old_cds = sum(1 for f in old.features if f.type == "CDS")
        new_cds = sum(1 for f in new.features if f.type == "CDS")
        print(f"{lname:<8}  {len(old.seq):>13,}  {len(new.seq):>13,}  "
              f"{old_cds:>8}  {new_cds:>8}  {acc}")
    print(f"\nKL303: unchanged (no discriminating representative found in NCBI)")
    print(f"\nNext: python scripts/type_normalized.py \\")
    print(f"        --db DB/EC-K-typing_all_groups_v0.5.gbk --suffix v0.5norm --threads 8")


if __name__ == "__main__":
    main()
