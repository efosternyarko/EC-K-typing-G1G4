# EC-K-typing Group 1 & 4

A reference database of *Escherichia coli* capsule (K-antigen) loci for **Group 1 and Group 4** capsule types.

This database complements the [EC-K-typing](https://github.com/rgladstone/EC-K-typing) Group 2 & 3 database by Rebecca Gladstone, enabling comprehensive capsule typing across all four *E. coli* capsule groups.

> **Pre-release status:** This database is under active development and uses a **0.x versioning scheme** until it reaches production quality. The current release is **v0.9**. Versions will be numbered 0.1, 0.2, 0.3, etc. until the database is used in a peer-reviewed publication, at which point it will be released as **v1.0**.
>
> **Normalised scoring:** A post-processing step re-ranks Kaptive's raw scores by `AS / total_expected_gene_bp`, eliminating bitscore accumulation bias from large reference loci. This achieves **651/651 (100%) self-typing** and **100% typeability** across all 1,126 validation assemblies.
>
> **⚠️ Do not use the combined all-groups database for typing.** Prior to v0.7, a combined G1-4 + G2-3 all-groups database (`EC-K-typing_all_groups_vX.X.gbk`) was distributed. Validation against independent cohorts revealed a *wzy*-interference artefact: G2/G3 isolates carry a *wzy*-dependent O-antigen locus (present in all *E. coli*, regardless of capsule group) that scores above the genuine *kps*-encoded capsule locus when Kaptive searches both databases simultaneously. This causes genuine G2/G3 isolates to appear untypeable when the combined database is used. **The G2/G3 and G1/G4 databases must always be run sequentially** (see [Usage](#usage) below). The all-groups files are retained in `DB/` for reference but should not be used for typing.

## Background

*E. coli* capsular polysaccharides are classified into four groups based on their biosynthetic pathways:

| Group | Pathway | Gene locus | Existing database |
|-------|---------|------------|-------------------|
| Group 2 | ABC transporter | near *serA* | [EC-K-typing](https://github.com/rgladstone/EC-K-typing) |
| Group 3 | ABC transporter | near *serA* | [EC-K-typing](https://github.com/rgladstone/EC-K-typing) |
| Group 1 | Wzy-dependent | *cps* locus (near *his*) | **This database** |
| Group 4 | Wzy-dependent | *cps* locus (near *his*) | **This database** |

Groups 1 and 4 share the Wzy-dependent polymerisation pathway. Their capsule biosynthesis genes are located at the *cps* locus, flanked by *galF* (upstream) and *gnd* (downstream), with the *wza-wzb-wzc* export genes further upstream.

## Usage

G1/G4 and G2/G3 capsule groups are mutually exclusive — no *E. coli* isolate carries both. The two databases must therefore be run **sequentially**, not simultaneously. Running both in a single Kaptive search causes a *wzy*-interference artefact (see warning above) that renders many genuine G2/G3 isolates untypeable.

### Recommended two-step workflow

**Step 1 — G2/G3 typing** using the Gladstone EC-K-typing database:

```bash
kaptive assembly -a assemblies/*.fasta \
    -k DB/EC-K-typing_group2and3_v3.0.0.gbk \
    -o results_G23/
```

**Step 2 — G1/G4 typing** on the untypeables from Step 1, using normalised scoring (recommended):

```bash
# Collect untypeable assembly paths from Step 1
# (assemblies where Kaptive reported None or Untypeable)
ls untypeable/*.fasta > untypeable_list.txt

# Run Kaptive in --scores mode against the G1/G4 database
kaptive assembly -a untypeable/*.fasta \
    -k DB/EC-K-typing_group1and4_v0.9.gbk \
    --scores results_G14/kaptive_scores_G14.tsv \
    -t 8
```

Combine Step 1 (typeable) and Step 2 results for a complete capsule-type assignment across all four groups.

### Alternative: kpsM pre-screen

For large collections where running the full G2/G3 Kaptive database is time-prohibitive, you can replace Step 1 with a faster `kpsM` gene screen. The `kpsM` gene encodes an inner membrane component of the Group 2/3 ABC-transporter capsule export system — it is present in all Group 2/3 strains and absent from Group 1/4 strains, making it a reliable marker to separate the two groups.

> **Limitation:** The kpsM screen identifies whether an assembly is Group 2/3 or a G1/G4 candidate, but does **not** assign a G2/G3 K-locus type. Use this approach only when G2/G3 K-locus typing is not needed.

**Requirements:** `blastn` and `makeblastdb` (NCBI BLAST+) on your `PATH`.

```bash
# Screen a directory of assemblies; run Kaptive on kpsM-negative ones:
python3 scripts/kpsm_screen.py \
    --assembly-dir assemblies/ \
    --kpsm-ref DB/kpsM_reference.fasta \
    --kaptive-db DB/EC-K-typing_group1and4_v0.9.gbk \
    --output-dir results/kpsm_screen/ \
    --run-kaptive
```

`DB/kpsM_reference.fasta` contains the reference sequences for `kpsM` (Group 2 marker) and `kpsM_3` (Group 3 marker). Detection of either at ≥90% identity and ≥80% query coverage marks an assembly as Group 2/3; all others are passed to G1/G4 typing.

Assemblies where `kpsM_screen_results.tsv` reports `group = Group2_3` are Group 2/3 strains — G1/G4 typing should not be applied to them. Assemblies with `group = G1_G4_candidate` should be passed to Step 2.

### Why normalised scoring?

Standard Kaptive accumulates raw alignment score across all reference genes. Because large loci accumulate proportionally more score, smaller loci of the same KX type are systematically outscored — even at equal per-base identity. Standard mode achieves only ~43% typeability on G1/G4 assemblies; normalised scoring achieves **100%**.

Normalised scoring divides the raw alignment score by the total expected CDS length of each reference locus, converting it to a per-base identity metric that is size-independent:

```
normalised score (AS_norm)  =  AS  /  total_expected_gene_bp
```

A value of ~2.0 indicates perfect 100% identity across all reference genes.

### Applying normalised scoring to your scores matrix

After running `kaptive assembly --scores`, apply normalisation with this Python script:

```python
#!/usr/bin/env python3
"""
normalise_kaptive_scores.py  —  re-rank Kaptive scores by AS / total_expected_gene_bp

Usage:
    python normalise_kaptive_scores.py \
        --db  DB/EC-K-typing_group1and4_v0.9.gbk \
        --in  results_G14/kaptive_scores_G14.tsv \
        --out results_G14/kaptive_results_norm.tsv
"""
import argparse, re
import pandas as pd
from Bio import SeqIO
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db",  required=True)
    ap.add_argument("--in",  required=True, dest="scores_tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-coverage", type=float, default=0.50,
                    help="Min gene coverage to call Typeable (default 0.50)")
    args = ap.parse_args()

    # Compute total CDS bp per locus from the GenBank
    bp = {}
    for rec in SeqIO.parse(args.db, "genbank"):
        lname = rec.name.split("_")[0]
        for f in rec.features:
            if f.type == "source":
                for note in f.qualifiers.get("note", []):
                    m = re.search(r"K locus:\s*(\S+)", note)
                    if m:
                        lname = m.group(1)
        cds = [f for f in rec.features if f.type == "CDS"]
        bp[lname] = sum(len(f) for f in cds)

    df = pd.read_csv(args.scores_tsv, sep="\t")
    asm_col = df.columns[0]

    df["total_bp"]  = df["Locus"].map(bp).fillna(df["q_len"])
    df["AS_norm"]   = df["AS"] / df["total_bp"].clip(lower=1)

    rows = []
    for asm, grp in df.groupby(asm_col):
        active = grp[grp["AS"] > 0]
        if active.empty:
            rows.append({"Assembly": asm, "Best match locus": "none",
                         "Best match confidence": "Untypeable", "Genes found": 0,
                         "Genes expected": 0, "Gene coverage": "0.0%",
                         "Raw AS": 0, "Norm AS": 0.0})
            continue
        best = active.sort_values(["AS_norm","AS"], ascending=[False,False]).iloc[0]
        gf   = int(best["genes_found"])
        ge   = int(best["genes_expected"])
        cov  = gf / ge if ge > 0 else 0
        conf = "Typeable" if cov >= args.min_coverage else "Untypeable"
        rows.append({"Assembly": asm, "Best match locus": best["Locus"],
                     "Best match confidence": conf, "Genes found": gf,
                     "Genes expected": ge, "Gene coverage": f"{100*cov:.1f}%",
                     "Raw AS": int(best["AS"]), "Norm AS": round(best["AS_norm"],4)})

    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print(f"Written: {args.out}  ({len(rows)} assemblies)")

if __name__ == "__main__":
    main()
```

Save this as `normalise_kaptive_scores.py` and run:

```bash
python normalise_kaptive_scores.py \
    --db DB/EC-K-typing_group1and4_v0.9.gbk \
    --in results_G14/kaptive_scores_G14.tsv \
    --out results_G14/kaptive_results_norm.tsv
```

**Example output (`kaptive_results_norm.tsv`):**

The `Assembly` column contains the filename stem as reported by Kaptive (path and extension stripped). The `Norm AS` column is what matters for interpreting typing results — values close to 2.0 indicate a high-confidence match.

```
Assembly            Best match locus  Best match confidence  Genes found  Genes expected  Gene coverage  Raw AS  Norm AS
sample_A            KL302             Typeable               38           40              95.0%          68207   1.8059
sample_B            KL304             Typeable               44           44              100.0%         93072   2.0
sample_C            KL306             Typeable               37           37              100.0%         75810   2.0
sample_D            KL436             Typeable               16           16              100.0%         29576   1.9545
sample_E            KL305             Typeable               50           50              100.0%         88182   2.0
```

> Note: `Raw AS` can range from ~10,000 to >100,000 depending on locus size. It is **not** a percentage. The column that is comparable across all loci is `Norm AS`.

**Interpreting Norm AS:**

| Norm AS | Interpretation |
|---------|---------------|
| ~2.00 | Perfect: 100% identity across all reference genes |
| ≥ 1.90 | High confidence: all or nearly all genes found at high identity |
| 1.50 – 1.90 | Moderate: partial locus or some gene-level divergence |
| < 1.50 | Low confidence: possibly novel, highly divergent, or mixed type |

**Output columns:**

| Column | Description |
|--------|-------------|
| `Best match locus` | KL locus designation (e.g., `KL302`, `K24`) |
| `Best match confidence` | `Typeable` if ≥ 50% of expected genes found; `Untypeable` otherwise |
| `Genes found` / `Genes expected` | Count of reference CDS found vs expected |
| `Gene coverage` | Fraction of expected genes found (found ÷ expected) |
| `Raw AS` | Cumulative alignment score — **not** comparable across loci of different sizes (range ~10,000–100,000+) |
| `Norm AS` | Raw AS ÷ total expected CDS base-pairs — size-independent per-base identity (~2.0 = perfect) |

---

### Understanding the two Kaptive modes

Kaptive can be run in two fundamentally different modes, which are **independent** of each other and produce separate outputs:

**Standard mode** (`kaptive assembly -a genome.fasta -k db.gbk -o out.tsv`):

1. minimap2 aligns the full reference locus sequence against the assembly → locate the locus *region*
2. Within that region, align each reference CDS → count genes found
3. If step 1 fails (no clear full-sequence hit), the locus is never located and the assembly is reported as **None** or **Untypeable**, regardless of whether individual genes exist somewhere in the assembly

Standard Kaptive reports a graded confidence level based on the fraction and identity of genes found:

| Confidence | Meaning |
|-----------|---------|
| Perfect | 100% expected genes found at high identity |
| Very High | Nearly all genes; very high identity |
| High | Most genes found; some minor divergence |
| Good | Majority of genes found |
| Low | Few genes found; locus may be incomplete or divergent |
| None | No significant match |
| Untypeable | Gene coverage below 50% threshold |

**Scores mode** (`kaptive assembly -a genome.fasta -k db.gbk --scores matrix.tsv`), used for normalised scoring:

1. minimap2 aligns each reference CDS **directly against the whole assembly** — no locus region is identified first
2. Scores are accumulated per locus across all genes found anywhere in the assembly
3. Every assembly receives a score for every reference locus; the best-scoring locus is the type call

Because scores mode does not depend on first locating a locus region, it succeeds even when the locus is fragmented across multiple contigs or too divergent for a full-sequence hit. This is why **normalised scoring achieves 100% typeability** compared to **~43% for standard Kaptive** on the same assemblies.

The two modes are run separately and do not interact. Standard Kaptive does not fall back to scores mode when step 1 fails, and normalised scoring does not use or override standard Kaptive results. They should be thought of as complementary tools: standard Kaptive confirms locus structure and location; normalised scoring provides the most sensitive and unbiased type assignment.

---

### With standard Kaptive only (no normalisation)

If you only need locus position information and Kaptive's confidence grades, and do not need 100% typeability:

```bash
kaptive assembly -a untypeable/*.fasta \
    -k DB/EC-K-typing_group1and4_v0.9.gbk \
    -o results_G14/
```

**Example output (first 8 columns):**

```
Assembly          Best match locus  Best match type  Match confidence  Problems  Identity  Coverage   Length discrepancy
sample1.fasta     KL302             KL302            Typeable          !         99.1%     93.1%      n/a
sample2.fasta     KL336             KL336            Typeable          +!        99.7%     100.0%     -325 bp
sample3.fasta     KL127             KL127            Typeable          -         99.9%     100.0%     +12 bp
```

`Best match locus` and `Best match type` are the same for G1/G4 loci — the KL designation is both the locus name and the type name, since most G1/G4 loci do not yet have a classical K-antigen designation. The `Problems` column uses codes from Kaptive (e.g., `!` = genes found outside the locus; `+` = extra unexpected genes; `-` = missing expected genes; `?` = other issues).

Standard Kaptive will correctly type most Group 2 & 3 loci without any post-processing. For Group 1 & 4 loci, the normalised scoring workflow above is strongly recommended, as ~57% of assemblies will be marked Untypeable by standard Kaptive due to locus fragmentation or size-bias in the scoring.

---

### Note on `scripts/type_normalized.py`

The repository includes `scripts/type_normalized.py`, which is the **internal validation script** used to generate the pre-computed validation results in `DB/`. This script uses hardcoded paths to the developer's genome collection (BSI assemblies, NCBI representative genomes, and ATB novel locus FASTAs) and is **not** intended for direct use by external users.

To type your own genomes, use `kaptive assembly --scores` followed by the `normalise_kaptive_scores.py` snippet above, which is a fully self-contained normalisation script requiring only the database GenBank file and the scores TSV.

---

## Database

### Versions

| Version | G1/G4 loci | Source genomes | Self-typing | Notes |
|---------|-----------|----------------|-------------|-------|
| v0.1 | 46 (KL300–KL343) | 222 (subset accessible) | 46/46 (100%) | Initial release |
| v0.2 | 125 (KL300–KL423) | 1,112 BSI | 70/93 (75.3%) | All BSI no-hit isolates |
| v0.3 | 93 filtered (positional names) | 1,112 BSI | 70/93 (75.3%) standard; 88/93 (94.6%) normalised | Positional gene naming + normalised scoring |
| v0.3.1 | 93 filtered | 1,112 BSI + 592 NNS | 592/592 (100%) NNS typeable | KL388/KL391 updated with longer NNS representatives |
| v0.4 | 93 filtered | 565 BSI | 89/93 (95.7%) normalised | Conserved CDS stripped (galF, gnd, ugd, wza, wzb, wzc) |
| v0.5 | 93 filtered | 565 BSI + 2 NCBI | 91/93 (97.8%) normalised | KL306 and KL307 replaced with NCBI representatives |
| v0.6 | 93 filtered | 565 BSI + 2 NCBI | 91/93 (97.8%) normalised | Conserved flanking/export CDS restored (KL301 retains stripped annotations) |
| v0.7 | 93 filtered | 565 BSI + 2 NCBI | 91/93 (97.8%) normalised | Combined all-groups database discontinued; sequential workflow mandatory |
| v0.8 | 667 (93 BSI + 574 novel ATB) | 565 BSI + 2 NCBI + 574 ATB | — | AllTheBacteria expansion: KL392–KL968 added (initial release, pre-validation) |
| **v0.9** | **651 (92 BSI + 559 novel ATB)** | **567 BSI/NCBI + 559 ATB = 1,126** | **651/651 (100%)** | **Removed 15 artefact/oversized loci, KL337 (no unique representative), and 3 indistinguishable pairs (merged); 100% self-typing and typeability** |

### Reference loci (v0.9)

The v0.9 database contains **651 reference K-loci**:

- **92 BSI original loci** (K24, KL300–KL391, excluding KL337): extracted from *E. coli* bloodstream infection isolates with no hit to the Group 2/3 database
  - Size range: 30.4 – 56.2 kb
  - KX types covered: KX01, KX17, KX28, KX31, KX34, KX36, KX63, KX67, KX72, K24
- **559 novel ATB loci** (KL392–KL968, minus removed/merged): representatives of novel K-locus clusters identified by screening ~2.4 M genomes from [AllTheBacteria](https://doi.org/10.1101/2024.03.08.584059)
  - Clustered at 95% nucleotide identity, 80% bidirectional coverage
  - Sorted by cluster size (largest first → lowest KL number)
  - Size range: ~10 – 60 kb (variable region only)

### Files

#### v0.9 (current)

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.9.gbk` | **G1/G4 database (651 loci) — use as the second step after G2/G3 typing** |
| `DB/novel_rep_kl_map.tsv` | KL → representative ATB assembly mapping (559 novel loci) |
| `DB/novel_kl_summary.tsv` | Novel loci summary: KL, representative assembly, sequence length, cluster size |
| `DB/novel_rep_fastas/` | Representative FASTA assemblies for novel ATB loci (downloaded via `download_novel_rep_fastas.py`) |
| `DB/kaptive_scores_v0.9norm_full.tsv` | Full locus × assembly score matrix (651 loci × 1,126 assemblies) |
| `DB/kaptive_validation_results_v0.9norm_full.tsv` | Normalised typing results (651 loci × 1,126 assemblies) |
| `DB/kaptive_validation_summary_v0.9norm_full.tsv` | Per-assembly comparison table (normalised scoring) |
| `DB/KL_G1G4_mapping.tsv` | KL nomenclature mapping for 92 BSI loci (KL name, KX origin, source assembly, length) |
| `DB/KL_G1G4_mapping_filtered.tsv` | Mapping for 92 BSI loci used in self-typing validation |
| `scripts/download_novel_rep_fastas.py` | Downloads 559 novel ATB representative assemblies (ATB S3 → ENA → NCBI fallback) |
| `scripts/remove_oversized_loci.py` | Removes novel loci with >60 CDS (flanking sequence extraction artefacts) |
| `scripts/type_normalized.py` | Internal validation script (hardcoded paths; see [Usage](#usage) above) |

#### v0.8

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.8.gbk` | G1/G4 database (667 loci): initial AllTheBacteria expansion (pre-validation cleanup) |

#### v0.7

> The G1/G4 GenBank file is unchanged from v0.6. This version removes the combined all-groups database from the recommended workflow following identification of a *wzy*-interference artefact (see warning above).

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.6.gbk` | G1/G4 database (93 loci) — superseded by v0.9 |
| `DB/EC-K-typing_all_groups_v0.6.gbk` | Combined all-groups database (183 loci) — **retained for reference only; do not use for typing** |
| `DB/kaptive_scores_v0.6norm.tsv` | Full locus × assembly score matrix (93 loci × 567 assemblies) |
| `DB/kaptive_validation_results_v0.6norm.tsv` | Normalised typing results (AS / total_expected_gene_bp) |
| `DB/kaptive_validation_summary_v0.6norm.tsv` | Per-assembly comparison table (normalised scoring) |

#### v0.6

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.6.gbk` | G1/G4 database: conserved flanking/export CDS restored to 92/93 loci (KL301 retains stripped annotations) |
| `DB/EC-K-typing_all_groups_v0.6.gbk` | Combined all-groups database (183 loci) — superseded; see v0.7 warning above |
| `DB/kaptive_scores_v0.6norm.tsv` | Full locus × assembly score matrix (183 loci × 567 assemblies) |
| `DB/kaptive_validation_results_v0.6norm.tsv` | Normalised typing results (AS / total_expected_gene_bp) |
| `DB/kaptive_validation_summary_v0.6norm.tsv` | Per-assembly comparison table (normalised scoring) |
| `scripts/make_v06_db.py` | Restores conserved CDS to all G1/G4 loci except KL301; uses blastn liftover for KL306/KL307 (new sequences) |

#### v0.5

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.5.gbk` | G1/G4 database: conserved CDS stripped + KL306/KL307 replaced with better NCBI representatives |
| `DB/EC-K-typing_all_groups_v0.5.gbk` | Combined all-groups database (183 loci) |
| `DB/kaptive_scores_v0.5norm.tsv` | Full locus × assembly score matrix (183 loci × 567 assemblies: 565 BSI + CP099041 + CP070103) |
| `DB/kaptive_validation_results_v0.5norm.tsv` | Normalised typing results (AS / total_expected_gene_bp) |
| `DB/kaptive_validation_summary_v0.5norm.tsv` | Per-assembly comparison table (normalised scoring) |
| `scripts/make_v05_db.py` | Script that replaces KL306/KL307 with NCBI candidates via BLAST-based liftover |
| `scripts/blast_ncbi_candidates.py` | NCBI megaBLAST search for candidate assemblies per locus |
| `scripts/test_blast_candidates.py` | Downloads and Kaptive-tests candidate assemblies |
| `scripts/add_ktype_notes.py` | Adds `K type:KLxxx` notes to G1/G4 GenBank source features so Kaptive reports type names instead of `unknown (KLxxx)` |

#### v0.4

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.4.gbk` | G1/G4 database with conserved CDS stripped (galF, galF_2, gnd, ugd, wza, wzb, wzc removed from all 93 loci) |
| `DB/EC-K-typing_all_groups_v0.4.gbk` | Combined all-groups database (183 loci) — G2/G3 unchanged, G1/G4 with variable region only |
| `DB/kaptive_validation_results_v0.4.tsv` | Standard Kaptive typing results (565 BSI assemblies, v0.4 DB) |
| `DB/kaptive_scores_v0.4norm.tsv` | Full locus × assembly score matrix (183 loci × 565 assemblies, v0.4 DB) |
| `DB/kaptive_validation_results_v0.4norm.tsv` | Normalised typing results (AS / total_expected_gene_bp) |
| `DB/kaptive_validation_summary_v0.4norm.tsv` | Per-assembly comparison table (normalised scoring) |
| `scripts/make_v04_db.py` | Script that strips conserved CDS from v0.3.1 G1/G4 loci to produce v0.4 |

#### v0.3.1

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.3.1.gbk` | **G1/G4 database with updated KL388 (+6.6 kb, Mlw) and KL391 (+10.5 kb, Barnards) representatives** |
| `DB/EC-K-typing_all_groups_v0.3.1.gbk` | **Combined all-groups database (183 loci) with updated reps — superseded; do not use** |
| `DB/nns_rep_update_summary.tsv` | Annotation summary for the two updated representatives |
| `DB/nns_kaptive_results_v0.3.1.tsv` | Standard Kaptive typing output (592 NNS assemblies, v0.3.1 DB) |
| `DB/nns_kaptive_scores_v0.3.1.tsv` | Full locus × assembly score matrix (183 loci × 592 assemblies, v0.3.1 DB) |
| `DB/nns_kaptive_results_v0.3.1_norm.tsv` | Normalised typing results with collection labels (v0.3.1 DB) |
| `DB/nns_kaptive_summary_v0.3.1_norm.tsv` | Per-collection summary table (v0.3.1 DB) |

#### v0.3

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.3.gbk` | 93 filtered loci with systematic positional gene names (GenBank, Kaptive-compatible) |
| `DB/EC-K-typing_all_groups_v0.3.gbk` | Combined database — all 4 capsule groups (183 loci, Kaptive-ready) |
| `DB/kaptive_validation_results_v0.3.tsv` | Kaptive v3.1.0 typing results for 565 v0.3 source assemblies (standard scoring) |
| `DB/kaptive_validation_summary_v0.3.tsv` | Per-assembly summary with expected vs observed KL (v0.3, standard scoring) |
| `DB/kaptive_scores_v0.3norm.tsv` | Full locus × assembly score matrix from `kaptive --scores` (183 loci × 565 assemblies) |
| `DB/kaptive_validation_results_v0.3norm.tsv` | Normalised typing results (AS / total\_expected\_gene\_bp) |
| `DB/kaptive_validation_summary_v0.3norm.tsv` | Per-assembly comparison table (normalised scoring) |
| `scripts/type_normalized.py` | **Internal validation script** — re-ranks Kaptive scores by `AS / total_expected_gene_bp` |

#### v0.2 (retained for reference)

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.2.fasta` | All 125 reference locus sequences (FASTA) |
| `DB/EC-K-typing_group1and4_v0.2.gbk` | Annotated reference loci (GenBank, Kaptive-compatible) |
| `DB/EC-K-typing_group1and4_v0.2_filtered.fasta` | Filtered set of 93 loci ≥ 30 kb (FASTA) |
| `DB/EC-K-typing_all_groups_v0.2.gbk` | Combined database — all 4 capsule groups (183 loci, Kaptive-ready) |
| `DB/KL_G1G4_mapping.tsv` | KL nomenclature mapping for all 92 BSI loci (KL name, KX origin, source assembly, length) |
| `DB/KL_G1G4_mapping_filtered.tsv` | Mapping for the 92 BSI locus filtered set |
| `DB/cluster_info.tsv` | Clustering details (cluster members, representative sequences) |
| `DB/kaptive_validation_results_v0.2.tsv` | Kaptive v3.1.0 typing results for 565 v0.2 source assemblies |
| `DB/kaptive_validation_summary_v0.2.tsv` | Per-assembly summary (v0.2) |

#### v0.1 (retained for reference)

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.1.fasta` | Reference locus sequences (FASTA) |
| `DB/EC-K-typing_group1and4_v0.1.gbk` | Annotated reference loci (GenBank, Kaptive-compatible) |
| `DB/EC-K-typing_all_groups_v0.1.gbk` | Combined database — all 4 capsule groups (136 loci) |
| `DB/kaptive_validation_results_v0.1.tsv` | Kaptive v3.1.0 typing results for 222 v0.1 source assemblies |

#### Other files

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group2and3_v3.0.0.gbk` | Gladstone Group 2 & 3 database (90 loci, included for convenience) |
| `scripts/build_G1G4_db.py` | Pipeline script for locus extraction and database construction |
| `scripts/annotate_loci.py` | Annotation script (pyrodigal + Klebsiella K-locus BLASTp) |
| `scripts/name_loci_positional.py` | Positional gene naming script (v0.3) |
| `scripts/type_normalized.py` | Internal validation script — re-ranks Kaptive scores by `AS / total_expected_gene_bp` (see Usage above) |
| `scripts/download_novel_rep_fastas.py` | Downloads representative FASTA assemblies for novel ATB loci |
| `scripts/remove_oversized_loci.py` | Removes novel loci with >60 CDS (flanking sequence artefacts) |
| `flanking_genes/flanking_genes.fasta` | Flanking gene marker sequences used for locus detection |

### Nomenclature

| KL designation | Description |
|----------------|-------------|
| K24 | Direct K-type call (known capsule type) |
| KL300–KL391 | BSI loci identified from EnteroBase (KL337 removed — no discriminating public representative) |
| KL392–KL968 | Novel loci identified from AllTheBacteria screening (559 loci retained after cleanup) |

## Methods

### Locus extraction pipeline

1. **Flanking gene detection:** BLAST search for *galF*, *gnd*, *ugd*, *wza*, and *wzc* from *E. coli* K-12 MG1655 against genome assemblies
2. **Locus boundary definition:** The *cps* region between *galF* and *gnd* (including *wza-wzb-wzc* upstream) with 500 bp flanking sequence
3. **Fragmented assembly handling:** For genomes where *galF* and *gnd* are on different contigs, sequences are extracted from both contigs and concatenated with an N-spacer
4. **Clustering:** All-vs-all BLAST at 95% identity and 80% query coverage, greedy clustering by sequence length
5. **Representative selection:** Longest sequence per cluster (v0.1); manually replaced with median-length alternatives for oversized loci in v0.2 to prevent scoring bias in Kaptive

### AllTheBacteria expansion (v0.8/v0.9)

The v0.8/v0.9 database incorporates novel K-loci identified by screening the complete [AllTheBacteria](https://doi.org/10.1101/2024.03.08.584059) (ATB) collection (~2.4 M prokaryotic assemblies) with the v0.7 G1/G4 database in scores mode. Assemblies with no hit to the G2/G3 database and a G1/G4 normalised score below the thresholds for known loci were flagged as potentially novel. Novel locus sequences were clustered by [MMseqs2](https://github.com/soedinglab/MMseqs2) at 95% nucleotide identity and 80% bidirectional coverage (cluster mode 0, greedy set cover), yielding 577 novel clusters. Representatives were assigned new KL numbers starting at KL392 (sorted largest cluster first), giving KL392–KL968.

Quality control for v0.9 removed:
- **15 artefact loci** (12 with >60 CDS, containing extensive flanking chromosomal sequence; plus KL486, KL562, KL812 with other structural issues)
- **KL337** (BSI locus with no discriminating public representative — both cluster members self-type to KL389 or KL767)
- **3 indistinguishable pairs** merged: KL443/KL446 → KL443; KL716/KL968 → KL716; KL830/KL843 → KL830 (all genes of both loci found at 100% identity in each other's assemblies)

### Source data

The BSI reference loci were derived from 6,673 *E. coli* bloodstream infection (BSI) isolates originally extracted from [EnteroBase](https://enterobase.warwick.ac.uk/):

1. **Group 2/3 screening:** All 6,673 BSI genomes were typed against the [Gladstone Group 2 & 3 database](https://github.com/rgladstone/EC-K-typing). **1,112 isolates** had no hit, indicating they likely carry Group 1 or Group 4 capsule loci.
2. **FastKaptive typing:** The 1,112 no-hit isolates were typed with [FastKaptive](https://github.com/rmostowy/fastKaptive), which assigned KX-type designations across 15 distinct types (KX36: 382, KX34: 329, KX17: 195, KX01: 73, KX67: 41, KX31: 29, KX72: 27, and others).
3. **Genome access:** All **1,112 genomes** were obtained for locus extraction in v0.2, resolving the EnteroBase accession mappings that limited v0.1 to 222 isolates.
4. **Locus extraction and clustering:** The extraction pipeline (see above) was applied to all 1,112 genomes, yielding **125 distinct K-locus clusters** (KL300–KL423, plus K24 and K96) covering 10 KX types. A filtered set of **93 loci ≥ 30 kb** was used in the combined all-groups database, excluding partial or fragmented extractions.

### Annotation

#### v0.2 annotation

Gene prediction and annotation follows a two-stage approach:

1. **Gene prediction:** [pyrodigal](https://github.com/althonos/pyrodigal) in metagenomic mode identifies open reading frames across each locus
2. **Gene name transfer:** Predicted proteins are searched against the [Kaptive *Klebsiella* K-locus reference database](https://github.com/klebgenomics/Kaptive) via BLASTp (>=30% identity, >=50% coverage). Matching genes receive functional names from the *Klebsiella* reference (e.g., *wza*, *wzb*, *wzc*, *wzx*, *wzy*, glycosyltransferases). 73% of predicted CDS received gene name annotations. Unnamed CDS receive their `locus_tag` as gene name (e.g., `KL302_00005`).

#### v0.3 positional naming

In v0.3, all predicted proteins across the 93 filtered loci are clustered by all-vs-all BLASTp (≥90% identity, ≥90% query and subject coverage), and family names are assigned deterministically:

- **Functional name:** if any cluster member has a Klebsiella-derived name (e.g., *wza*, *gmd*), the most common such name is propagated to the whole family
- **Positional name:** otherwise the family is named `KL{N}_{pos}`, where N is the KL number of the lowest-numbered locus that contains a member, and pos is its ordinal position in that locus — identical across every locus that shares the protein family

Result: 3,298 CDS total; 2,094 with functional names (64%); 1,204 with positional names (36%); 67 multi-locus protein families with a shared name. Within-locus duplicate names receive `_2`, `_3` suffixes.

GenBank records are formatted for [Kaptive](https://github.com/klebgenomics/Kaptive) compatibility, with each locus identified by a `source` feature containing `/note="K locus: KLxxx"`.

## Validation

### v0.9 validation

The v0.9 database was validated by running `kaptive assembly --scores` on all 1,126 representative assemblies (567 BSI/NCBI + 559 ATB novel) against the 651-locus database and applying normalised scoring:

| Metric | Result |
|--------|--------|
| Self-typing (651 loci, all representatives) | **651/651 (100.0%)** |
| Self-typing — 92 original BSI loci | **92/92 (100.0%)** |
| Self-typing — 559 novel ATB loci | **559/559 (100.0%)** |
| Assemblies typeable | **1,126/1,126 (100.0%)** |

Full validation results: `DB/kaptive_validation_results_v0.9norm_full.tsv` · `DB/kaptive_validation_summary_v0.9norm_full.tsv`

### v0.1 validation

The v0.1 combined database was validated by re-typing the 222 v0.1 source genome assemblies with Kaptive v3.1.0:

| Metric | Result |
|--------|--------|
| Self-typing (46 reference loci) | **46/46 (100%)** |
| Genomes typeable | **161/222 (72.5%)** |
| All 46 G1/4 loci utilised | Yes |
| Untypeable genomes | 61 (mostly Group 2/3 types without G1/4 locus) |

Full validation results: `DB/kaptive_validation_results_v0.1.tsv`

### v0.2 validation

The v0.2 database was validated by running Kaptive v3.1.0 on all 565 assemblies from which a locus was successfully extracted:

| Metric | Result |
|--------|--------|
| Self-typing (93 filtered loci) | **70/93 (75.3%)** |
| Self-typing (125 full set) | 70/125 (56.0%) |
| Assemblies typeable | **258/565 (45.7%)** |
| Loci utilised (filtered set) | 70/93 |

The 23 loci that do not self-type are predominantly KX01 loci that score below KL302 (40 CDS, 41 kb). See v0.3 validation for root-cause analysis.

Full validation results: `DB/kaptive_validation_results_v0.2.tsv`

### v0.3 validation

The v0.3 database (positional gene naming applied) was validated identically to v0.2:

| Metric | Result |
|--------|--------|
| Self-typing (93 filtered loci) | **70/93 (75.3%)** |
| Self-typing (125 full set) | 70/125 (56.0%) |
| Assemblies typeable | 258/565 (45.7%) |
| Loci utilised (filtered set) | 70/93 |

Self-typing is **unchanged** from v0.2. Positional gene naming did not improve Kaptive's discrimination because Kaptive scores by cumulative BLAST bitscore across all reference genes — gene names in the GenBank qualifiers do not affect scoring. The root cause of the 23 failures is **bitscore accumulation bias from oversized representatives**: KL302 has 40 CDS (37.8 kb total CDS length) and KL304 has 44 CDS (46.5 kb), accumulating more total bitscore than smaller within-KX01/KX34 references (28–35 kb) even when all query genes match at 100% identity.

Full validation results: `DB/kaptive_validation_results_v0.3.tsv`

### Normalised scoring (v0.3norm)

`scripts/type_normalized.py` post-processes the full Kaptive scores matrix by re-ranking each assembly's loci by `AS / total_expected_gene_bp` (alignment score per expected reference base). This converts raw bitscore accumulation into a per-base identity metric that is size-independent:

| Metric | Standard Kaptive | Normalised scoring |
|--------|------------------|--------------------|
| Self-typing (93 filtered loci) | 70/93 (75.3%) | **88/93 (94.6%)** |
| Self-typing (125 full set) | 70/125 (56.0%) | 88/125 (70.4%) |
| Assemblies typeable | 258/565 (45.7%) | **565/565 (100.0%)** |
| Loci utilised (filtered set) | 70/93 | 88/93 |

**Normalisation fixes 20 of the 23 failures** (KL300, KL301, KL303, KL306, KL307 still fail — all KX01 loci that genuinely score lower per base than KL302 in their source genomes, likely reflecting biological similarity within this KX type).

Full results: `DB/kaptive_validation_results_v0.3norm.tsv` · `DB/kaptive_validation_summary_v0.3norm.tsv`

### With BLAST

For direct BLAST-based typing:

```bash
makeblastdb -in DB/EC-K-typing_group1and4_v0.2_filtered.fasta -dbtype nucl

blastn -query genome.fasta \
  -db DB/EC-K-typing_group1and4_v0.2_filtered.fasta \
  -outfmt '6 qseqid sseqid pident qcovs length' \
  -max_target_seqs 5
```

## Neonatal sepsis (NNS) typing

The database was applied to **592 *E. coli* assemblies** from 8 neonatal sepsis collections representing 6 countries/sites, using the v0.5 all-groups database with normalised scoring:

| Collection | Site | n | G1/G4 (%) | Top G1/G4 loci |
|------------|------|:-:|:---------:|----------------|
| Patan | Nepal | 20 | 75% (15) | KL302(8), KL329(5), KL337(1) |
| Barnards | South Africa | 75 | 64% (48) | KL302(24), KL301(5), KL303(5) |
| Mbira | Zimbabwe | 57 | 46% (26) | KL303(6), KL329(5), KL302(4) |
| CHAMPS_Harar | Ethiopia | 110 | 54% (59) | KL302(19), KL301(13), KL303(8) |
| Mlw | Malawi | 167 | 48% (80) | KL302(32), KL301(13), KL303(5) |
| MRCG | Gambia | 130 | 34% (44) | KL302(15), KL379(9), KL305(5) |
| Benin | Benin | 20 | 55% (11) | KL329(4), KL302(4), KL303(1) |
| Pakistan | Pakistan | 13 | 62% (8) | KL302(4), KL303(2), KL344(1) |

**All 592/592 assemblies were typeable** (v0.5 database, normalised scoring). G1/G4 loci were assigned to 291 assemblies (49.2%); G2/G3 loci to 301 assemblies (50.8%). The most common G1/G4 type was KL302 (110/291; 38%), followed by KL301 (35), KL303 (28), KL329 (23), and KL337 (19). Full per-assembly results: `DB/nns_kaptive_results_v0.5norm.tsv`.

### Database improvement candidates from NNS

40 G1/G4 loci were represented in the NNS collections. High-quality matches (normalised AS ≥ 1.9, 100% gene coverage) were found for 22 loci including KL305, KL309, KL313, KL315, KL317, KL319, KL323, KL324, KL327, KL329, KL332, KL336, KL344, KL361, KL363, KL365, KL371, KL374, KL379, KL381, KL388, and KL391. The best NNS representative per locus is catalogued in `DB/nns_g1g4_best_reps.tsv`. KL388 and KL391 were updated in v0.3.1.

### NNS typing files (`DB/`)

| File | Description |
|------|-------------|
| `nns_kaptive_results_v0.5norm.tsv` | **Normalised typing results, v0.5 database** — Assembly, Collection, Best_match_locus, Capsule_group (G1/G4 or G2/G3), Confidence (High/Moderate/Low), Norm_AS, Genes_found, Genes_expected, Gene_coverage, Raw_AS, Typeable |
| `nns_kaptive_scores_v0.5norm.tsv` | Full locus × assembly score matrix (183 loci × 592 assemblies, v0.5 DB) |
| `nns_kaptive_results.tsv` | Standard Kaptive typing output (592 assemblies, 8 collections, v0.3 DB) |
| `nns_kaptive_results_norm.tsv` | Normalised typing results with collection labels (v0.3 DB) |
| `nns_kaptive_summary_norm.tsv` | Per-collection summary table (v0.3 DB) |
| `nns_kaptive_scores.tsv` | Full locus × assembly score matrix (183 loci × 592 assemblies, v0.3 DB) |
| `nns_g1g4_best_reps.tsv` | Best NNS assembly per G1/G4 locus (DB expansion candidates) |

## Roadmap

### v0.3 (completed)

Implemented systematic positional gene naming across all 93 filtered loci and score normalisation post-processing:

**Positional naming** (`scripts/name_loci_positional.py`):
- **Conserved genes:** retain functional Klebsiella-derived names (*wza*, *wzb*, *wzc*, *wzx*, *wzy*, *galF*, *gnd*, *ugd*, glycosyltransferases)
- **Variable-region genes:** clustered at 90% identity/coverage across all loci; shared families receive a single positional name (`KL{N}_{pos}`) named from the lowest-numbered locus that contributes a member
- **Result:** 3,298 CDS; 64% functional names; 36% positional names; 67 multi-locus shared protein families

**Finding:** positional naming did not change standard Kaptive self-typing (70/93, 75.3%). Gene qualifier names in GenBank records do not affect Kaptive's raw bitscore accumulation.

**Score normalisation** (`scripts/type_normalized.py`):
- Re-ranks Kaptive's per-locus alignment scores by `AS / total_expected_gene_bp` — alignment score per expected reference base
- Eliminates bitscore accumulation bias from large reference loci (KL302: 37.8 kb CDS; KL304: 46.5 kb CDS)
- **Result: 88/93 (94.6%) — fixes 20 of 23 failures; typeability 100%**
- 5 remaining failures: KL300, KL301, KL303, KL306, KL307 → all KX01 loci that genuinely score lower per base than KL302

### v0.3.1 (completed)

Replaced two representatives with longer, more complete sequences extracted from neonatal sepsis (NNS) assemblies:

- **KL388:** 33,951 bp (BSI) → 40,565 bp (+6.6 kb; Malawi/Mlw, ERR11525392); 36 CDS, 100% annotated
- **KL391:** 30,360 bp (BSI) → 40,871 bp (+10.5 kb; South Africa/Barnards, ERR4920077); 36 CDS, 97% annotated

**Key fix:** KL391 assemblies previously mistyped as KL300 (Untypeable) now correctly type as KL391 (Typeable) with standard Kaptive scoring.

**NNS re-validation (v0.3.1):** All 592/592 NNS assemblies remain typeable with normalised scoring. ERR4920077 and ERR4920086 (Barnards, South Africa) now type as KL391 with 100% gene coverage (NormAS 2.00 and 1.88 respectively) — both previously typed as KL391 via normalised scoring on v0.3 and now also correctly typed by standard Kaptive. NNS validation files: `DB/nns_kaptive_results_v0.3.1_norm.tsv`, `DB/nns_kaptive_summary_v0.3.1_norm.tsv`.

### v0.4 (completed)

Stripped conserved flanking/export CDS from all 93 G1/G4 GenBank records, leaving only the variable biosynthetic region for Kaptive scoring. The 5 KX01 loci that failed in v0.3.1 (KL300, KL301, KL303, KL306, KL307) all shared identical conserved genes (wza, wzb, wzc, galF, gnd, ugd) with KL302, contributing equal bitscore background that masked variable region differences.

**Implementation** (`scripts/make_v04_db.py`):
- Read v0.3.1 G1/G4 GenBank; remove any CDS whose `gene` qualifier is in `{galF, galF_2, gnd, ugd, wza, wzb, wzc}`
- G2/G3 records are carried through unchanged
- 528 / 3,304 CDS (16.0%) removed across 93 loci; each locus retains 23–47 variable CDS (above Kaptive's 50% gene coverage threshold)

**Result:**
- Self-typing (93 filtered loci): **89/93 (95.7%)** — up from 88/93 (94.6%)
- KL301 now correctly self-types (previously failed)
- Typeability: **565/565 (100%)** — unchanged
- Remaining 4 failures: KL300, KL303, KL306, KL307 (all KX01; still score below KL302 in variable region)

**Interpretation:** The conserved gene removal partially helps — KL301 is resolved — but the 4 remaining failures reflect genuine similarity of the KX01 variable biosynthetic genes, not scoring bias from conserved background. Further improvement likely requires better representative sequences or additional discriminating genes.

### v0.5 (completed)

Replaced the representative sequences for KL300, KL306, and KL307 with better candidates identified by NCBI megaBLAST and validated with Kaptive normalised scoring. KL303 was left unchanged — no publicly available *E. coli* genome was found whose KL303 locus discriminates from KL302 with normalised scoring.

**Candidate discovery** (`scripts/blast_ncbi_candidates.py`):
- Submitted each failing locus sequence to NCBI megaBLAST (nt, *Escherichia coli* taxid:562, ≥90% identity, ≥70% query coverage)
- KL300: 299 hits | KL303: 264 hits | KL306: 6 hits | KL307: 30 hits

**Candidate validation** (`scripts/test_blast_candidates.py`):
- Downloaded chromosome FASTAs for the top 10 candidates per locus via NCBI Entrez
- Ran Kaptive `--scores` against the v0.4 database and applied normalised scoring
- KL300: 10/10 correct | KL303: 0/10 (all → KL302 or KL352) | KL306: 6/6 correct | KL307: 4/10 correct

**Replacement representatives** (`scripts/make_v05_db.py`):
- Locus extracted from chromosome by multi-HSP BLAST spanning; annotation transferred by per-CDS blastn liftover (no re-annotation with pyrodigal)
- Conserved CDS stripped as in v0.4

| Locus | Old rep (BSI) | New rep | Old len | New len | CDS mapped |
|-------|--------------|---------|---------|---------|------------|
| KL306 | ESC_GB3606AA_AS | CP099041 | 41,735 bp | 42,626 bp | 31/31 |
| KL307 | ESC_NB5901AA_AS | CP070103 | 44,329 bp | 45,329 bp | 34/34 |

**KL300** (CP135488) was evaluated but not adopted. With the 1.2× extraction cap (27/30 CDS), CP135488's core variable genes are too conserved across KX01 loci — the new reference pulled many non-KL300 assemblies toward KL300 in normalised scoring. With a wider 1.5× cap (all 30 CDS including peripheral rfaG/gmd_2/KL300_11 flanking genes), standard Kaptive typeability regressed by 54 assemblies due to size-bias in locus-finding. KL300 retains its original BSI representative.

**Self-typing validation methodology correction:** Self-typing should use the source genome of each locus's *current* representative. KL306's representative was extracted from CP099041 and KL307's from CP070103; those genomes are therefore the correct self-typing assemblies. The mapping tables (`KL_G1G4_mapping.tsv`, `KL_G1G4_mapping_filtered.tsv`) were updated accordingly.

**Result:**
- Self-typing (93 filtered loci, correct source genomes): **91/93 (97.8%) normalised** — up from 89/93 with legacy BSI test set
- Correctly self-typing: all 93 loci except **KL300** (old BSI rep types as KL302) and **KL303** (biological ambiguity with KL302)
- New NCBI representatives self-type at **AS_norm = 2.000** (perfect): CP099041 → KL306, CP070103 → KL307
- Typeability (normalised): **567/567 (100%)** — includes 565 BSI + 2 NCBI assemblies

**KL303 note:** All 264 NCBI candidate assemblies typed as KL302 or KL352 — the KL303 variable region is indistinguishable from KL302 in all publicly available sequences. KL303 is retained in the database but flagged as having no discriminating public representative.

### v0.6 (completed)

Restored conserved flanking/export CDS annotations (wza, wzb, wzc, galF, galF_2, gnd, ugd) to all G1/G4 GenBank records except KL301. These genes were stripped in v0.4 to reduce scoring bias, but validation showed the stripping fixed only KL301 — the other failures (KL306, KL307) were resolved by representative replacement in v0.5, while KL300 and KL303 remain biologically unresolvable regardless of whether conserved genes are present.

**Rationale for restoration:** Retaining full locus annotations (including conserved flanking genes) is important for correct gene-level characterisation of assembled loci — partial annotations could mislead users about what genes are present in a K-locus. Normalised scoring performance is identical with or without conserved genes for all loci except KL301.

**Implementation** (`scripts/make_v06_db.py`):
- **KL301:** conserved CDS kept stripped — this is what resolved it in v0.4
- **91 loci with unchanged sequences (same as v0.3.1):** conserved CDS features copied directly from the v0.3.1 GenBank — coordinates are identical
- **KL306 and KL307 (new NCBI representatives):** conserved CDS positions located in the new sequences by per-gene blastn (≥80% identity, ≥80% query coverage), with length snapped to a multiple of 3

**Result:**
- Self-typing (93 filtered loci): **91/93 (97.8%) normalised** — identical to v0.5
- Typeability: **567/567 (100%)**
- Remaining failures: **KL300** and **KL303** (KX01 loci biologically indistinguishable from KL302)

### v0.7 (completed)

Combined all-groups database (`EC-K-typing_all_groups_vX.X.gbk`) discontinued after validation against Norwegian, UK, and Malawi cohorts demonstrated a *wzy*-interference artefact causing G2/G3 isolates to be mis-scored as untypeable when G1/G4 and G2/G3 databases are searched simultaneously. Sequential two-step workflow (G2/G3 first, then G1/G4 on untypeables) is now the mandatory approach. G1/G4 database unchanged from v0.6.

### v0.8 (completed)

AllTheBacteria (ATB) expansion: screened ~2.4 M prokaryotic assemblies with the v0.7 database in scores mode, identified assemblies not matching any known K-locus at high confidence, extracted novel K-locus sequences, and clustered with MMseqs2 (95% identity, 80% bidirectional coverage). Assigned new KL numbers KL392–KL968 (577 novel clusters, sorted largest-first). Initial release with 93 BSI + 574 novel ATB = **667 loci** total.

Remaining failure from v0.7: **KL300** and **KL303** (KX01 loci with no discriminating public representative). KL300 and KL303 were accepted as known ambiguous loci in the context of the full ATB expansion.

### v0.9 (completed)

Quality control pass on v0.8 to remove loci causing systematic typing failures:

1. **Oversized loci removed** (12 novel loci with >60 CDS): extraction pipeline captured extensive flanking chromosomal sequence; these scored ~2.0 in every assembly, overriding the correct type assignment
2. **Three additional problematic novel loci removed** (KL486, KL562, KL812): structural issues identified during validation
3. **KL337 removed**: both cluster members (only 2 total) self-type to KL389 or KL767, not KL337; no assembly in the 1,126-genome validation set types uniquely to KL337
4. **Three indistinguishable pairs merged**: KL443/KL446 → KL443; KL716/KL968 → KL716; KL830/KL843 → KL830 (all genes of both loci found at 100% identity in each other's assemblies)
5. **Tiebreaker improved**: replaced old 0.5% window + genes_found tiebreaker with pure AS_norm descending then raw AS descending, eliminating systematic bias toward larger loci in near-tie cases

**Result: 651/651 (100%) self-typing; 1,126/1,126 (100%) typeability**

### v1.0 (public release)

v1.0 will be released when:
- The database has been used in a peer-reviewed publication
- Systematic gene naming is implemented for novel ATB loci (currently using locus_tag placeholders)
- AllTheBacteria expansion is updated with the latest ATB release

## Related work

| Tool/Database | Scope | Reference |
|---------------|-------|-----------|
| [EC-K-typing](https://github.com/rgladstone/EC-K-typing) | *E. coli* Group 2 & 3 (90 loci) | [Gladstone et al. 2024](https://www.medrxiv.org/content/10.1101/2024.11.22.24317484v1) |
| [kTYPr](https://github.com/SushiLab/kTYPr) | *E. coli* Group 2 & 3 (85 loci, HMM-based) | [Schwengers et al. 2025](https://www.biorxiv.org/content/10.1101/2025.08.07.669119v1) |
| [Kaptive](https://github.com/klebgenomics/Kaptive) | *Klebsiella* and *E. coli* K/O typing | [Lam et al. 2022](https://doi.org/10.1099/mgen.0.000800) |
| [FastKaptive](https://github.com/rmostowy/fastKaptive) | Fast K-locus pre-screening | Mostowy et al. |
| **This database** | *E. coli* Group 1 & 4 (651 loci, v0.9; 651/651 with normalised scoring) | — |

## Citation

If you use this database, please cite:

- This repository
- Gladstone RA et al. (2024). *E. coli* capsule typing. [medRxiv preprint](https://www.medrxiv.org/content/10.1101/2024.11.22.24317484v1)
- The [EC-K-typing](https://github.com/rgladstone/EC-K-typing) Group 2 & 3 database

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
