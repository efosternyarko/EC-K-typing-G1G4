# EC-K-typing Group 1 & 4

A reference database of *Escherichia coli* capsule (K-antigen) loci for **Group 1 and Group 4** capsule types.

This database complements the [EC-K-typing](https://github.com/rgladstone/EC-K-typing) Group 2 & 3 database by Rebecca Gladstone, enabling comprehensive capsule typing across all four *E. coli* capsule groups.

> **Pre-release status:** This database is under active development and uses a **0.x versioning scheme** until it reaches production quality. The current release is **v0.5**. Versions will be numbered 0.1, 0.2, 0.3, etc. until the database passes full self-typing validation and systematic gene naming is implemented, at which point it will be released as **v1.0**.
>
> **Normalised scoring:** A post-processing script (`scripts/type_normalized.py`) re-ranks Kaptive's raw scores by `AS / total_expected_gene_bp`, eliminating bitscore accumulation bias from large reference loci. This raises self-typing of the 93 filtered loci from **70/93 (75.3%)** to **91/93 (97.8%)** and typeability to **100%** of isolates.

## Background

*E. coli* capsular polysaccharides are classified into four groups based on their biosynthetic pathways:

| Group | Pathway | Gene locus | Existing database |
|-------|---------|------------|-------------------|
| Group 2 | ABC transporter | near *serA* | [EC-K-typing](https://github.com/rgladstone/EC-K-typing) |
| Group 3 | ABC transporter | near *serA* | [EC-K-typing](https://github.com/rgladstone/EC-K-typing) |
| Group 1 | Wzy-dependent | *cps* locus (near *his*) | **This database** |
| Group 4 | Wzy-dependent | *cps* locus (near *his*) | **This database** |

Groups 1 and 4 share the Wzy-dependent polymerisation pathway. Their capsule biosynthesis genes are located at the *cps* locus, flanked by *galF* (upstream) and *gnd* (downstream), with the *wza-wzb-wzc* export genes further upstream.

## Database

### Versions

| Version | G1/G4 loci | All-groups loci | Source genomes | Self-typing | Notes |
|---------|-----------|-----------------|----------------|-------------|-------|
| v0.1 | 46 (KL300–KL343) | 136 | 222 (subset accessible) | 46/46 (100%) | Initial release |
| v0.2 | 125 (KL300–KL423) | 183 | 1,112 (all BSI no-hit isolates) | 70/93 (75.3%) | |
| v0.3 | 93 filtered (positional names) | 183 | 1,112 BSI | 70/93 (75.3%) standard; 88/93 (94.6%) normalised | Positional gene naming + normalised scoring script |
| v0.3.1 | 93 filtered | 183 | 1,112 BSI + 592 NNS | 592/592 (100%) NNS typeable | KL388 (+6.6 kb) and KL391 (+10.5 kb) replaced with longer NNS representatives from Malawi and South Africa; KL391 now correctly typed by standard Kaptive scoring (previously KL300 Untypeable) |
| v0.4 | 93 filtered | 183 | 565 BSI | 89/93 (95.7%) normalised | Conserved CDS stripped from G1/G4 loci (528 CDS removed: galF, galF_2, gnd, ugd, wza, wzb, wzc); fixes KL301; typeability 100% |
| **v0.5** | **93 filtered** | **183** | **565 BSI + 2 NCBI** | **91/93 (97.8%) normalised; 67/93 (72.0%) standard** | **KL306 and KL307 representatives replaced with better NCBI candidates (CP099041, CP070103); self-typing now uses source genomes of current representatives; new reps self-type at AS_norm=2.00; typeability 100%** |

### Reference loci (v0.2)

The v0.2 database contains **125 reference K-loci** (K24, K96, KL300–KL423) extracted from all 1,112 *E. coli* bloodstream infection isolates with no hit to the Group 2/3 database:

- **Size range (full set):** 7.0 – 59.1 kb
- **Total sequence (full set):** 4.14 Mb
- **Filtered set (≥ 30 kb):** 93 loci, 3.69 Mb — used in the combined all-groups database
- **Clustering threshold:** 95% nucleotide identity, 80% query coverage
- **KX types covered:** KX01, KX17, KX28, KX31, KX34, KX36, KX63, KX67, KX72, K24

### Files

#### v0.5 (current)

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.5.gbk` | **G1/G4 database: conserved CDS stripped + KL300/KL306/KL307 replaced with better NCBI representatives** |
| `DB/EC-K-typing_all_groups_v0.5.gbk` | **Combined all-groups database (183 loci) — use this for typing** |
| `DB/kaptive_scores_v0.5norm.tsv` | Full locus × assembly score matrix (183 loci × 567 assemblies: 565 BSI + CP099041 + CP070103) |
| `DB/kaptive_validation_results_v0.5norm.tsv` | Normalised typing results (AS / total_expected_gene_bp) |
| `DB/kaptive_validation_summary_v0.5norm.tsv` | Per-assembly comparison table (normalised scoring) |
| `scripts/make_v05_db.py` | Script that replaces KL306/KL307 with NCBI candidates via BLAST-based liftover |
| `scripts/blast_ncbi_candidates.py` | NCBI megaBLAST search for candidate assemblies per locus |
| `scripts/test_blast_candidates.py` | Downloads and Kaptive-tests candidate assemblies |

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
| `DB/EC-K-typing_all_groups_v0.3.1.gbk` | **Combined all-groups database (183 loci) with updated reps — use this for typing** |
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
| `scripts/type_normalized.py` | **Normalised typing script** — re-ranks Kaptive scores by `AS / total_expected_gene_bp` |

#### v0.2 (retained for reference)

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v0.2.fasta` | All 125 reference locus sequences (FASTA) |
| `DB/EC-K-typing_group1and4_v0.2.gbk` | Annotated reference loci (GenBank, Kaptive-compatible) |
| `DB/EC-K-typing_group1and4_v0.2_filtered.fasta` | Filtered set of 93 loci ≥ 30 kb (FASTA) |
| `DB/EC-K-typing_all_groups_v0.2.gbk` | Combined database — all 4 capsule groups (183 loci, Kaptive-ready) |
| `DB/KL_G1G4_mapping.tsv` | KL nomenclature mapping for all 125 loci (KL name, KX origin, source assembly, length) |
| `DB/KL_G1G4_mapping_filtered.tsv` | Mapping for the 93-locus filtered set |
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
| `scripts/type_normalized.py` | Normalised scoring script — re-ranks Kaptive scores by `AS / total_expected_gene_bp` |
| `scripts/swap_reps_v04.py` | Representative swapping script (v0.4 experiment — kept for reference, superseded by make_v04_db.py) |
| `scripts/make_v04_db.py` | v0.4 database builder — strips conserved CDS from G1/G4 loci for variable-region-only scoring |
| `flanking_genes/flanking_genes.fasta` | Flanking gene marker sequences used for locus detection |

### Nomenclature

| KL designation | Description |
|----------------|-------------|
| K24, K96 | Direct K-type calls (known capsule types) |
| KL300–KL423 | Novel locus types identified in this study |

## Methods

### Locus extraction pipeline

1. **Flanking gene detection:** BLAST search for *galF*, *gnd*, *ugd*, *wza*, and *wzc* from *E. coli* K-12 MG1655 against genome assemblies
2. **Locus boundary definition:** The *cps* region between *galF* and *gnd* (including *wza-wzb-wzc* upstream) with 500 bp flanking sequence
3. **Fragmented assembly handling:** For genomes where *galF* and *gnd* are on different contigs, sequences are extracted from both contigs and concatenated with an N-spacer
4. **Clustering:** All-vs-all BLAST at 95% identity and 80% query coverage, greedy clustering by sequence length
5. **Representative selection:** Longest sequence per cluster (v0.1); manually replaced with median-length alternatives for oversized loci in v0.2 to prevent scoring bias in Kaptive

### Source data

The reference loci were derived from 6,673 *E. coli* bloodstream infection (BSI) isolates originally extracted from [EnteroBase](https://enterobase.warwick.ac.uk/):

1. **Group 2/3 screening:** All 6,673 BSI genomes were typed against the [Gladstone Group 2 & 3 database](https://github.com/rgladstone/EC-K-typing). **1,112 isolates** had no hit, indicating they likely carry Group 1 or Group 4 capsule loci.
2. **FastKaptive typing:** The 1,112 no-hit isolates were typed with [FastKaptive](https://github.com/rmostowy/fastKaptive), which assigned KX-type designations across 15 distinct types (KX36: 382, KX34: 329, KX17: 195, KX01: 73, KX67: 41, KX31: 29, KX72: 27, and others).
3. **Genome access:** All **1,112 genomes** were obtained for locus extraction in v0.2, resolving the EnteroBase accession mappings that limited v0.1 to 222 isolates.
4. **Locus extraction and clustering:** The extraction pipeline (see above) was applied to all 1,112 genomes, yielding **125 distinct K-locus clusters** (KL300–KL423, plus K24 and K96) covering 10 KX types. A filtered set of **93 loci ≥ 30 kb** is provided for use with Kaptive, excluding partial or fragmented extractions.

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

## Usage

### With Kaptive

The current combined database (`DB/EC-K-typing_all_groups_v0.5.gbk`) covering all four capsule groups (183 loci) is ready for direct use with [Kaptive](https://github.com/klebgenomics/Kaptive):

```bash
kaptive assembly DB/EC-K-typing_all_groups_v0.5.gbk genome.fasta -o results.tsv
```

The combined database merges:
- **Group 2 & 3:** 90 loci (KL1–KL175) from [Gladstone et al.](https://github.com/rgladstone/EC-K-typing) — annotated with Bakta + Panaroo
- **Group 1 & 4:** 93 loci (K24, K96, KL300–KL423, filtered ≥ 30 kb) from this study — annotated with [pyrodigal](https://github.com/althonos/pyrodigal) (metagenomic mode) + systematic positional gene naming (v0.3); KL388 and KL391 updated to longer NNS representatives (v0.3.1); conserved flanking/export genes stripped for variable-region-only scoring (v0.4); KL306 and KL307 replaced with better NCBI representatives (v0.5)

### With normalised scoring (recommended)

For improved accuracy on Group 1 & 4 loci, use the normalised scoring script, which re-ranks each assembly's loci by `AS / total_expected_gene_bp` — converting raw bitscore into alignment score per expected reference base:

```bash
# Run Kaptive and normalise scores in one step (v0.5 database):
python scripts/type_normalized.py \
  --db DB/EC-K-typing_all_groups_v0.5.gbk \
  --suffix v0.5norm \
  --threads 8

# Or re-use a pre-computed scores matrix (fast):
python scripts/type_normalized.py --skip-kaptive --suffix v0.5norm
```

This raises self-typing of the 93 filtered loci to **91/93 (97.8%)** and typeability to **100%** of assemblies. See `DB/kaptive_validation_results_v0.5norm.tsv` for the pre-computed results on all 567 assemblies (565 BSI + 2 NCBI representative genomes).

## Validation

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

The near-tie tiebreaker — among loci scoring within 0.5% of the best normalised score, prefer the locus with more genes found — resolves one additional ambiguous case (KL337 vs KL389, 0.17% AS_norm difference, 41 vs 29 genes found).

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

The database was applied to **592 *E. coli* assemblies** from 8 neonatal sepsis collections representing 6 countries/sites:

| Collection | Site | n assemblies | G1/G4 (%) | Top G1/G4 loci |
|------------|------|:------------:|:---------:|----------------|
| Patan | Nepal | 20 | 100% | KL302(16), KL329(10), KL337(2) |
| Barnards | South Africa | 75 | 40% | KL302(28), KL301(5), KL303(5) |
| Mbira | Zimbabwe | 57 | 42% | KL303(6), KL127(6), KL329(5) |
| CHAMPS_Harar | Ethiopia | 110 | 50% | KL302(17), KL301(11), KL303(7) |
| Mlw | Malawi | 167 | 43% | KL302(41), KL127(16), KL130(14) |
| MRCG | Gambia | 130 | 48% | KL2(28), KL130(17), KL302(12) |
| Benin | Benin | 20 | 55% | KL329(4), KL302(4), KL127(2) |
| Pakistan | Pakistan | 13 | 54% | KL302(4), KL303(2), KL131(2) |

**All 592/592 assemblies were typeable** using the all-groups v0.3 database with normalised scoring. G1/G4 loci were assigned to 317 assemblies (53.5%); G2/G3 loci to 308 assemblies (47.0%). The most common G1/G4 type was KL302 (126/317; 40%), followed by KL303 (28), KL329 (27), KL301 (21), and KL337 (20).

### Database improvement candidates from NNS

40 G1/G4 loci were represented in the NNS collections. High-quality matches (normalised AS ≥ 1.9, 100% gene coverage) were found for 22 loci including KL305, KL309, KL313, KL315, KL317, KL319, KL323, KL324, KL327, KL329, KL332, KL336, KL344, KL361, KL363, KL365, KL371, KL374, KL379, KL381, KL388, and KL391. The best NNS representative per locus is catalogued in `DB/nns_g1g4_best_reps.tsv`. KL388 and KL391 were updated in v0.3.1.

### NNS typing files (`DB/`)

| File | Description |
|------|-------------|
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
- Near-tie tiebreaker: among loci within 0.5% of top normalised score, prefer the one with more genes found
- **Result: 88/93 (94.6%) — fixes 20 of 23 failures; typeability 100%**
- 5 remaining failures: KL300, KL301, KL303, KL306, KL307 → all KX01 loci that genuinely score lower per base than KL302

**Note on representative swapping (tested, not released):** Swapping KL302 and KL305 to smaller representatives was attempted (see `scripts/swap_reps_v04.py`) but yielded 64/93 (68.8%) — worse than v0.3. Bitscore accumulation bias is systemic; shrinking one dominant locus shifts dominance to the next-largest. Score normalisation is the correct fix.

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

**Self-typing validation methodology correction:** Self-typing should use the source genome of each locus's *current* representative. KL306's representative was extracted from CP099041 and KL307's from CP070103; those genomes are therefore the correct self-typing assemblies. The mapping tables (`KL_G1G4_mapping.tsv`, `KL_G1G4_mapping_filtered.tsv`) were updated accordingly, and `scripts/type_normalized.py` now also searches the `DB/blast_ncbi_results/candidate_genomes/` directory for representative genomes not in the original BSI assembly set. The old BSI assemblies (ESC_GB3606AA_AS, ESC_NB5901AA_AS) remain in the 565 BSI test set but are no longer designated as the KL306/KL307 representatives; they type as KL302 even with the v0.5 database, confirming that replacement was necessary.

**Result:**
- Self-typing (93 filtered loci, correct source genomes): **91/93 (97.8%) normalised** — up from 89/93 with legacy BSI test set
- Correctly self-typing: all 93 loci except **KL300** (old BSI rep types as KL302) and **KL303** (biological ambiguity with KL302)
- New NCBI representatives self-type at **AS_norm = 2.000** (perfect): CP099041 → KL306, CP070103 → KL307
- Standard Kaptive: **67/93 (72.0%) self-typing; 243/565 (43.0%) typeable** — identical to v0.4
- Typeability (normalised): **567/567 (100%)** — includes 565 BSI + 2 NCBI assemblies

**KL303 note:** All 264 NCBI candidate assemblies typed as KL302 or KL352 — the KL303 variable region is indistinguishable from KL302 in all publicly available sequences. KL303 is retained in the database but flagged as having no discriminating public representative.

### v0.6 (next release)

Remaining failures are **KL300** and **KL303** (both KX01 loci). KL303 and KL300 cannot be resolved from KL302 using any public NCBI genome. Options:

#### Option A: AllTheBacteria search

- Mine [AllTheBacteria](https://doi.org/10.1101/2024.03.08.584059) (~2.4M prokaryotic assemblies) using [LexicMap](https://www.nature.com/articles/s41587-025-02812-8) for discriminating representatives — requires ~64 GB RAM (EC2 c7g.8xlarge)

#### Option B: Accept biological ambiguity

- Document KL300, KL303, KL306, KL307 as known ambiguous loci; flag assemblies assigned to these types with a low-confidence warning in the output

### v1.0 (public release)

v1.0 will be released when:
- Self-typing reaches **100%** for all loci in the filtered set
- AllTheBacteria expansion is complete
- The database has been used in a peer-reviewed publication

## Related work

| Tool/Database | Scope | Reference |
|---------------|-------|-----------|
| [EC-K-typing](https://github.com/rgladstone/EC-K-typing) | *E. coli* Group 2 & 3 (90 loci) | [Gladstone et al. 2024](https://www.medrxiv.org/content/10.1101/2024.11.22.24317484v1) |
| [kTYPr](https://github.com/SushiLab/kTYPr) | *E. coli* Group 2 & 3 (85 loci, HMM-based) | [Schwengers et al. 2025](https://www.biorxiv.org/content/10.1101/2025.08.07.669119v1) |
| [Kaptive](https://github.com/klebgenomics/Kaptive) | *Klebsiella* and *E. coli* K/O typing | [Lam et al. 2022](https://doi.org/10.1099/mgen.0.000800) |
| [FastKaptive](https://github.com/rmostowy/fastKaptive) | Fast K-locus pre-screening | Mostowy et al. |
| **This database** | *E. coli* Group 1 & 4 (93 filtered loci, v0.5; 91/93 with normalised scoring) | — |

## Citation

If you use this database, please cite:

- This repository
- Gladstone RA et al. (2024). *E. coli* capsule typing. [medRxiv preprint](https://www.medrxiv.org/content/10.1101/2024.11.22.24317484v1)
- The [EC-K-typing](https://github.com/rgladstone/EC-K-typing) Group 2 & 3 database

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
