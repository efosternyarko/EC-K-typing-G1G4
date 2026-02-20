#!/usr/bin/env python3
"""
make_v06_db.py — Restore conserved CDS annotations to all G1/G4 loci except KL301.

Background
----------
v0.4 stripped wza/wzb/wzc/galF/galF_2/gnd/ugd from all 93 G1/G4 records to
remove scoring bias from shared flanking/export genes.  Validation showed that
removing them fixed only one locus (KL301), while the other four failures
(KL306, KL307: fixed by better NCBI representatives in v0.5; KL300, KL303:
biologically unresolvable) were fixed or remained unfixable regardless of
whether conserved genes are present.

Strategy
---------
  - KL301: keep conserved genes stripped (their removal is what fixed it).
  - 91 loci with unchanged sequences (same as v0.3.1): copy conserved CDS
    features directly from the v0.3.1 GenBank — coordinates are identical.
  - KL306 and KL307: new representative sequences (from CP099041 / CP070103);
    use per-gene blastn to lift conserved CDS positions from the v0.3.1
    KL306/KL307 records onto the new sequences.

Output
------
  DB/EC-K-typing_group1and4_v0.6.gbk   — G1/G4 only
  DB/EC-K-typing_all_groups_v0.6.gbk   — combined (G2/G3 unchanged + v0.6 G1/G4)
"""
import copy
import re
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

# ---------------------------------------------------------------------------
REPO_DIR  = Path(__file__).resolve().parent.parent
DB_DIR    = REPO_DIR / "DB"
V031_G1G4 = DB_DIR / "EC-K-typing_group1and4_v0.3.1.gbk"
V05_G1G4  = DB_DIR / "EC-K-typing_group1and4_v0.5.gbk"
G2G3_GBK  = DB_DIR / "EC-K-typing_group2and3_v3.0.0.gbk"
V06_G1G4  = DB_DIR / "EC-K-typing_group1and4_v0.6.gbk"
V06_ALL   = DB_DIR / "EC-K-typing_all_groups_v0.6.gbk"

CONSERVED_NAMES = {'galF', 'galF_2', 'gnd', 'ugd', 'wza', 'wzb', 'wzc'}
SKIP_LOCI       = {'KL301'}   # keep conserved genes stripped — that's the fix
LIFTOVER_LOCI   = {'KL306', 'KL307'}  # new sequences: use blastn

MIN_BLAST_PIDENT  = 80.0   # min % identity for blastn liftover
MIN_BLAST_QCOV    = 80.0   # min % query coverage — reject partial CDS
# ---------------------------------------------------------------------------


def get_locus_name(rec) -> str:
    for f in rec.features:
        if f.type == "source":
            for note in f.qualifiers.get("note", []):
                m = re.search(r"K locus:\s*(\S+)", note)
                if m:
                    return m.group(1)
    return rec.name.split("_")[0]


def conserved_features(rec) -> list:
    """Return conserved CDS features from a record."""
    return [
        f for f in rec.features
        if f.type == "CDS"
        and f.qualifiers.get("gene", [""])[0] in CONSERVED_NAMES
    ]


def blastn_lift(query_seq: str, subject_seq: str, gene_name: str,
                source_feat: SeqFeature) -> SeqFeature | None:
    """
    Align query_seq (one CDS from v0.3.1) against subject_seq (new rep).
    Returns a new SeqFeature at the found position, or None if no hit.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        q = Path(tmpdir) / "q.fa"
        s = Path(tmpdir) / "s.fa"
        q.write_text(f">q\n{query_seq}\n")
        s.write_text(f">s\n{subject_seq}\n")
        r = subprocess.run(
            ["blastn", "-query", str(q), "-subject", str(s),
             "-outfmt", "6 sstart send pident qcovs",
             "-perc_identity", str(MIN_BLAST_PIDENT),
             "-max_target_seqs", "1", "-max_hsps", "1"],
            capture_output=True, text=True
        )
    if not r.stdout.strip():
        return None
    parts = r.stdout.strip().split("\n")[0].split("\t")
    sstart, send, pident, qcovhsp = int(parts[0]), int(parts[1]), float(parts[2]), float(parts[3])
    if qcovhsp < MIN_BLAST_QCOV:
        print(f"      [{gene_name}] qcov={qcovhsp:.0f}% < {MIN_BLAST_QCOV:.0f}% — likely truncated at locus edge, skipping")
        return None
    if sstart <= send:
        start0, end0, strand = sstart - 1, send, 1
    else:
        start0, end0, strand = send - 1, sstart, -1
    # Snap to multiple of 3 (blastn can be off by ±1-2 bp at truncated boundaries)
    length = end0 - start0
    rem = length % 3
    if rem != 0:
        adj = 3 - rem  # 1 or 2 bp to add
        # prefer extending toward the interior of the sequence
        if strand == -1:   # gene reads rightward on minus strand: extend start
            start0 = max(0, start0 - adj)
        else:              # gene reads leftward on plus strand: extend end
            end0 = min(end0 + adj, len(subject_seq))
        length = end0 - start0
        if length % 3 != 0:   # still off — trim instead
            end0 -= length % 3
    loc = FeatureLocation(start0, end0, strand=strand)
    new_feat = SeqFeature(loc, type="CDS",
                          qualifiers=copy.deepcopy(source_feat.qualifiers))
    return new_feat


# ---------------------------------------------------------------------------
# 1. Load v0.3.1 conserved features indexed by locus name
# ---------------------------------------------------------------------------
print("[1] Loading v0.3.1 conserved CDS features...")
v031_conserved: dict[str, list] = {}
v031_recs: dict[str, object] = {}
for rec in SeqIO.parse(str(V031_G1G4), "genbank"):
    ln = get_locus_name(rec)
    feats = conserved_features(rec)
    v031_conserved[ln] = feats
    v031_recs[ln] = rec
    if feats:
        print(f"    {ln:8s}: {len(feats)} conserved CDS  ({[f.qualifiers.get('gene',[''])[0] for f in feats]})")

# ---------------------------------------------------------------------------
# 2. Process v0.5 G1/G4 records
# ---------------------------------------------------------------------------
print("\n[2] Restoring conserved CDS features...")
v05_recs = list(SeqIO.parse(str(V05_G1G4), "genbank"))
stats = {"skipped": 0, "direct_copy": 0, "blastn_lift": 0, "no_v031": 0}

for rec in v05_recs:
    ln = get_locus_name(rec)

    if ln in SKIP_LOCI:
        print(f"    {ln:8s}: SKIPPED (conserved genes stripped intentionally)")
        stats["skipped"] += 1
        continue

    if ln in LIFTOVER_LOCI:
        # New sequence: blast each conserved CDS from v0.3.1 onto new seq
        src = v031_recs.get(ln)
        if src is None:
            print(f"    {ln:8s}: WARNING — no v0.3.1 record found, skipping")
            stats["no_v031"] += 1
            continue
        added = []
        for feat in v031_conserved[ln]:
            gene = feat.qualifiers.get("gene", ["?"])[0]
            cds_seq = str(feat.extract(src.seq))
            new_feat = blastn_lift(cds_seq, str(rec.seq), gene, feat)
            if new_feat:
                added.append(new_feat)
                print(f"    {ln:8s}: {gene} → blastn hit at {new_feat.location}")
            else:
                print(f"    {ln:8s}: {gene} → NO BLASTN HIT (check manually)")
        rec.features.extend(added)
        stats["blastn_lift"] += 1

    else:
        # Same sequence: copy features directly (coordinates are identical)
        feats = v031_conserved.get(ln, [])
        if not feats:
            stats["no_v031"] += 1
            continue
        rec.features.extend([copy.deepcopy(f) for f in feats])
        stats["direct_copy"] += 1

    # Sort all features by position
    rec.features.sort(key=lambda f: (int(f.location.start), int(f.location.end)))

print(f"\n    Skipped (KL301): {stats['skipped']}")
print(f"    Direct copy (same seq): {stats['direct_copy']}")
print(f"    Blastn liftover (new seq): {stats['blastn_lift']}")
print(f"    No v0.3.1 record: {stats['no_v031']}")

# ---------------------------------------------------------------------------
# 3. Verify — print CDS counts before/after
# ---------------------------------------------------------------------------
print("\n[3] Verification (conserved CDS count in v0.6):")
for rec in v05_recs:
    ln = get_locus_name(rec)
    n_cons = len(conserved_features(rec))
    n_var  = len([f for f in rec.features
                  if f.type == "CDS" and
                  f.qualifiers.get("gene", [""])[0] not in CONSERVED_NAMES])
    n_expected = v031_conserved.get(ln, [])
    if ln in SKIP_LOCI:
        ok = "OK" if n_cons == 0 else "MISMATCH"
    elif ln in LIFTOVER_LOCI:
        ok = "OK (liftover)" if n_cons <= len(n_expected) else "MISMATCH"
    else:
        ok = "OK" if n_cons == len(n_expected) else "MISMATCH"
    print(f"    {ln:8s}: cons={n_cons}  var={n_var}  expected={len(n_expected)}  {ok}")

# ---------------------------------------------------------------------------
# 4. Write outputs
# ---------------------------------------------------------------------------
print("\n[4] Writing databases...")
SeqIO.write(v05_recs, str(V06_G1G4), "genbank")
print(f"    Written: {V06_G1G4.name}")

g2g3_recs = list(SeqIO.parse(str(G2G3_GBK), "genbank"))
SeqIO.write(v05_recs + g2g3_recs, str(V06_ALL), "genbank")
print(f"    Written: {V06_ALL.name}  ({len(v05_recs)} G1/G4 + {len(g2g3_recs)} G2/G3 = {len(v05_recs)+len(g2g3_recs)} total)")
print("\nDone. Next: run type_normalized.py --db DB/EC-K-typing_all_groups_v0.6.gbk --suffix v0.6norm")
