# Database Construction Protocol

How the G1/G4 K-locus reference database was built, from raw assemblies to the current versioned release.

---

## v1.2 — Discovery and annotation (623 loci)

**1. Identify G1/G4 candidates**
Bloodstream infection (BSI) *E. coli* assemblies from EnteroBase were typed against the Gladstone G2/G3 database using Kaptive. Assemblies that were untypeable (no G2/G3 locus) are G1/G4 candidates.

**2. Extract K-loci**
The K-locus region was extracted from each candidate assembly by identifying flanking conserved genes (*galF* upstream, *gnd* downstream). This yielded one locus sequence per assembly.

**3. Cluster and select representatives**
Loci were clustered at 80% nucleotide identity (≥95% coverage) using CD-HIT and MMseqs2. One representative was selected per cluster — the longest complete sequence within each cluster. This removed near-identical duplicates and retained one example of each distinct locus type.

**4. Annotate representatives**
Representative loci were annotated using Prokka and manually curated against known capsule biosynthesis gene names (wza, wzb, wzc, wzy, wzx, gtr*, rml*, etc.). Annotated records were compiled into a Kaptive-compatible GenBank file → **v1.2** (623 loci, KL300–KL967).

---

## v1.3 — Variant consolidation (490 loci)

v1.2 still contained loci that encode the same capsule type but differ by mobile element insertions, small gene deletions, or tandem duplications. These were identified and collapsed.

**5. Detect IS elements**
All locus proteins were screened against ISfinder HMM profiles ([`scripts/01_extract_proteins.py`](scripts/01_extract_proteins.py), [`scripts/02_detect_is_proteins.sh`](scripts/02_detect_is_proteins.sh)). IS-encoding CDS were flagged.

**6. Compare gene content with Panaroo**
All 623 loci were converted to GFF3 format ([`scripts/02b_gbk_to_gff3.py`](scripts/02b_gbk_to_gff3.py)) and analysed with Panaroo at 70% protein identity ([`scripts/02c_run_panaroo.sh`](scripts/02c_run_panaroo.sh)) to generate a pan-genome presence/absence matrix. IS gene families were excluded from the matrix ([`scripts/04_build_presence_absence.py`](scripts/04_build_presence_absence.py)).

**7. Identify variants**
Using the IS-free presence/absence matrix, four classes of variant were identified ([`scripts/04b_panaroo_find_variants.py`](scripts/04b_panaroo_find_variants.py)):

| Class | Definition | n removed |
|---|---|---|
| IS variant | Same biosynthetic gene content; one member has an IS element | 73 |
| Biosynthetic deletion | One locus is a strict gene-content subset of another (≤5 genes missing, zero unique genes) | 76 |
| Direct repeat variant | Intra-locus tandem duplication of a biosynthetic gene block (100% internal identity) | 1 |
| Small locus subset | Core biosynthetic genes are a subset of a larger parent; flanking gene differences only | 5 |

For IS variants, the IS-free locus is the primary; for all other classes, the gene-richer locus is the primary. Clinker visualisations of IS-variant groups are in [`DB/clinker_is_variants/`](DB/clinker_is_variants/).

**8. Build v1.3**
Variant loci were removed and IS CDS annotations were stripped from retained loci ([`scripts/06_build_v1.3.py`](scripts/06_build_v1.3.py)). The complete merge table is [`DB/is_analysis/panaroo_loci_to_merge.tsv`](DB/is_analysis/panaroo_loci_to_merge.tsv).

**9. Validate**
Each of the 490 loci in v1.3 was typed back against the database using Kaptive ([`scripts/07_kaptive_validation.py`](scripts/07_kaptive_validation.py)). All 490 self-type correctly (100%). Results: [`DB/kaptive_validation/`](DB/kaptive_validation/).

---

## Version summary

| Version | Loci | Notes |
|---|---|---|
| v1.2 | 623 | Discovery set; includes IS variants, deletion variants, direct repeats |
| v1.3 | 490 | Variant-consolidated; 100% kaptive self-typing |
