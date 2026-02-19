# EC-K-typing Group 1 & 4

A reference database of *Escherichia coli* capsule (K-antigen) loci for **Group 1 and Group 4** capsule types.

This database complements the [EC-K-typing](https://github.com/rgladstone/EC-K-typing) Group 2 & 3 database by Rebecca Gladstone, enabling comprehensive capsule typing across all four *E. coli* capsule groups.

> **Pre-release status:** This database is under active development and uses a **0.x versioning scheme** until it reaches production quality. The current release is **v0.3**. Versions will be numbered 0.1, 0.2, 0.3, etc. until the database passes full self-typing validation and systematic gene naming is implemented, at which point it will be released as **v1.0**. File names within the repository retain internal version numbers (v1.0, v2.0, v3.0) and will be reconciled at the v1.0 release.

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
| **v0.3** | **93 filtered (positional names)** | **183** | **1,112 BSI** | **70/93 (75.3%)** | **Current release — systematic positional gene naming implemented** |
| v0.4 | TBD | TBD | 1,112 BSI | Target: ≥90/93 | Planned — representative swapping for oversized loci |

### Reference loci (v2.0)

The v2.0 database contains **125 reference K-loci** (K24, K96, KL300–KL423) extracted from all 1,112 *E. coli* bloodstream infection isolates with no hit to the Group 2/3 database:

- **Size range (full set):** 7.0 – 59.1 kb
- **Total sequence (full set):** 4.14 Mb
- **Filtered set (≥ 30 kb):** 93 loci, 3.69 Mb — used in the combined all-groups database
- **Clustering threshold:** 95% nucleotide identity, 80% query coverage
- **KX types covered:** KX01, KX17, KX28, KX31, KX34, KX36, KX63, KX67, KX72, K24

### Files

#### v3.0 (current)

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v3.0.gbk` | **93 filtered loci with systematic positional gene names (GenBank, Kaptive-compatible)** |
| `DB/EC-K-typing_all_groups_v3.0.gbk` | **Combined database** — all 4 capsule groups (183 loci, Kaptive-ready) |
| `DB/kaptive_validation_results_v3.tsv` | Kaptive v3.1.0 typing results for 565 v3.0 source assemblies |
| `DB/kaptive_validation_summary_v3.tsv` | Per-assembly summary with expected vs observed KL (v3.0) |

#### v2.0 (retained for reference)

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v2.0.fasta` | All 125 reference locus sequences (FASTA) |
| `DB/EC-K-typing_group1and4_v2.0.gbk` | Annotated reference loci (GenBank, Kaptive-compatible) |
| `DB/EC-K-typing_group1and4_v2.0_filtered.fasta` | Filtered set of 93 loci ≥ 30 kb (FASTA) |
| `DB/EC-K-typing_all_groups_v2.0.gbk` | Combined database — all 4 capsule groups (183 loci, Kaptive-ready) |
| `DB/KL_G1G4_mapping.tsv` | KL nomenclature mapping for all 125 loci (KL name, KX origin, source assembly, length) |
| `DB/KL_G1G4_mapping_filtered.tsv` | Mapping for the 93-locus filtered set |
| `DB/cluster_info.tsv` | Clustering details (cluster members, representative sequences) |
| `DB/kaptive_validation_results_v2.tsv` | Kaptive v3.1.0 typing results for 565 v2.0 source assemblies |

#### v1.0 (retained for reference)

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v1.0.fasta` | Reference locus sequences (FASTA) |
| `DB/EC-K-typing_group1and4_v1.0.gbk` | Annotated reference loci (GenBank, Kaptive-compatible) |
| `DB/EC-K-typing_all_groups_v1.0.gbk` | Combined database — all 4 capsule groups (136 loci) |
| `DB/kaptive_validation_results.tsv` | Kaptive v3.1.0 typing results for 222 v1.0 source assemblies |

#### Other files

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group2and3_v3.0.0.gbk` | Gladstone Group 2 & 3 database (90 loci, included for convenience) |
| `scripts/build_G1G4_db.py` | Pipeline script for locus extraction and database construction |
| `scripts/annotate_loci.py` | Annotation script (pyrodigal + Klebsiella K-locus BLASTp) |
| `scripts/name_loci_positional.py` | Positional gene naming script (v3.0) |
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
3. **Genome access:** All **1,112 genomes** were obtained for locus extraction in v2.0, resolving the EnteroBase accession mappings that limited v1.0 to 222 isolates.
4. **Locus extraction and clustering:** The extraction pipeline (see above) was applied to all 1,112 genomes, yielding **125 distinct K-locus clusters** (KL300–KL423, plus K24 and K96) covering 10 KX types. A filtered set of **93 loci ≥ 30 kb** is provided for use with Kaptive, excluding partial or fragmented extractions.

### Annotation

#### v2.0 annotation

Gene prediction and annotation follows a two-stage approach:

1. **Gene prediction:** [pyrodigal](https://github.com/althonos/pyrodigal) in metagenomic mode identifies open reading frames across each locus
2. **Gene name transfer:** Predicted proteins are searched against the [Kaptive *Klebsiella* K-locus reference database](https://github.com/klebgenomics/Kaptive) via BLASTp (>=30% identity, >=50% coverage). Matching genes receive functional names from the *Klebsiella* reference (e.g., *wza*, *wzb*, *wzc*, *wzx*, *wzy*, glycosyltransferases). 73% of predicted CDS received gene name annotations. Unnamed CDS receive their `locus_tag` as gene name (e.g., `KL302_00005`).

#### v3.0 positional naming

In v3.0, all predicted proteins across the 93 filtered loci are clustered by all-vs-all BLASTp (≥90% identity, ≥90% query and subject coverage), and family names are assigned deterministically:

- **Functional name:** if any cluster member has a Klebsiella-derived name (e.g., *wza*, *gmd*), the most common such name is propagated to the whole family
- **Positional name:** otherwise the family is named `KL{N}_{pos}`, where N is the KL number of the lowest-numbered locus that contains a member, and pos is its ordinal position in that locus — identical across every locus that shares the protein family

Result: 3,298 CDS total; 2,094 with functional names (64%); 1,204 with positional names (36%); 67 multi-locus protein families with a shared name. Within-locus duplicate names receive `_2`, `_3` suffixes.

GenBank records are formatted for [Kaptive](https://github.com/klebgenomics/Kaptive) compatibility, with each locus identified by a `source` feature containing `/note="K locus: KLxxx"`.

## Usage

### With Kaptive

A pre-built combined database (`DB/EC-K-typing_all_groups_v3.0.gbk`) covering all four capsule groups (183 loci) is included and ready for direct use with [Kaptive](https://github.com/klebgenomics/Kaptive):

```bash
kaptive assembly DB/EC-K-typing_all_groups_v3.0.gbk genome.fasta -o results.tsv
```

The combined database merges:
- **Group 2 & 3:** 90 loci (KL1–KL175) from [Gladstone et al.](https://github.com/rgladstone/EC-K-typing) — annotated with Bakta + Panaroo
- **Group 1 & 4:** 93 loci (K24, K96, KL300–KL423, filtered ≥ 30 kb) from this study — annotated with [pyrodigal](https://github.com/althonos/pyrodigal) (metagenomic mode) + systematic positional gene naming (v3.0)

## Validation

### v1.0 validation

The v1.0 combined database was validated by re-typing the 222 v1.0 source genome assemblies with Kaptive v3.1.0:

| Metric | Result |
|--------|--------|
| Self-typing (46 reference loci) | **46/46 (100%)** |
| Genomes typeable | **161/222 (72.5%)** |
| All 46 G1/4 loci utilised | Yes |
| Untypeable genomes | 61 (mostly Group 2/3 types without G1/4 locus) |

Full validation results: `DB/kaptive_validation_results.tsv`

### v3.0 validation

The v3.0 database (positional gene naming applied) was validated identically to v2.0:

| Metric | Result |
|--------|--------|
| Self-typing (93 filtered loci) | **70/93 (75.3%)** |
| Self-typing (125 full set) | 70/125 (56.0%) |
| Assemblies typeable | 258/565 (45.7%) |
| Loci utilised (filtered set) | 70/93 |

Self-typing is **unchanged** from v2.0. Positional gene naming did not improve Kaptive's discrimination because Kaptive scores by cumulative BLAST bitscore across all reference genes — gene names in the GenBank qualifiers do not affect scoring. The root cause of the 23 failures is **bitscore accumulation bias from oversized representatives**: KL302 has 40 CDS (41 kb) and KL305 has 50 CDS (51 kb); these loci accumulate more total bitscore than smaller within-KX01 references (30–38 CDS) even when all query genes match at 100% identity. Fixing this requires swapping oversized representatives to size-matched alternatives (planned for v0.4).

Full validation results: `DB/kaptive_validation_results_v3.tsv`

### v2.0 validation

The v2.0 database was validated by running Kaptive v3.1.0 on all 565 assemblies from which a locus was successfully extracted:

| Metric | Result |
|--------|--------|
| Self-typing (93 filtered loci) | **70/93 (75.3%)** |
| Self-typing (125 full set) | 70/125 (56.0%) |
| Assemblies typeable | **258/565 (45.7%)** |
| Loci utilised (filtered set) | 70/93 |

The 23 loci that do not self-type are predominantly KX01 loci that score below KL302 (40 CDS, 41 kb). See v3.0 validation for root-cause analysis.

Full validation results: `DB/kaptive_validation_results_v2.tsv`

### With BLAST

For direct BLAST-based typing:

```bash
makeblastdb -in DB/EC-K-typing_group1and4_v2.0_filtered.fasta -dbtype nucl

blastn -query genome.fasta \
  -db DB/EC-K-typing_group1and4_v2.0_filtered.fasta \
  -outfmt '6 qseqid sseqid pident qcovs length' \
  -max_target_seqs 5
```

## Roadmap

### v0.3 (completed)

Implemented systematic positional gene naming across all 93 filtered loci (see `scripts/name_loci_positional.py`):

- **Conserved genes:** retain functional Klebsiella-derived names (*wza*, *wzb*, *wzc*, *wzx*, *wzy*, *galF*, *gnd*, *ugd*, glycosyltransferases)
- **Variable-region genes:** clustered at 90% identity/coverage across all loci; shared families receive a single positional name (`KL{N}_{pos}`) named from the lowest-numbered locus that contributes a member; singleton families get a name tied to their locus
- **Result:** 3,298 CDS; 64% functional names; 36% positional names; 67 multi-locus shared protein families

**Finding:** self-typing unchanged at 70/93 (75.3%). Positional naming does not affect Kaptive's scoring because Kaptive accumulates raw BLAST bitscore — gene qualifier names in GenBank records play no role. The 23 persistent failures are caused by bitscore accumulation bias: KL302 (40 CDS) and KL305 (50 CDS) are oversized representatives that out-score same-KX-type loci with 30–38 CDS regardless of per-gene identity.

### v0.4 (attempted — not released)

Representative swapping was attempted for KL302 and KL305, the two most dominant loci:

| Locus | Old representative | Old CDS | New representative | New CDS |
|-------|-------------------|---------|-------------------|---------|
| KL302 | ESC_GB6443AA_AS | 40 | ESC_CC1376AA_AS (KX01) | 28 |
| KL305 | ESC_BA8240AA_AS | 50 | ESC_TA7291AA_AS (KX34) | 37 |

**Result: 64/93 (68.8%) — worse than v3.0.** The swap fixed 6 loci that were being dominated by KL302 (KL326, KL337, KL364, KL365, KL367, KL368) but broke 12 others that were correctly self-typing in v3.0 (KL302 itself typed as KL309; many others shifted to being dominated by KL300 or KL304 instead). Net: −6 loci.

**Conclusion — whack-a-mole:** bitscore accumulation bias is systemic. Every KX type has a "most CDS" reference that accumulates the highest score in its neighbours' genomes. Shrinking one dominant locus shifts dominance to the next-largest. The 23 failures cannot be fixed one locus at a time.

### v0.4 (next release)

True resolution requires addressing the systemic scoring problem. Two viable approaches:

#### Option A: Score normalisation (preferred)

Add a Kaptive post-processing step that normalises each locus score by the number of expected genes (or locus length), converting raw bitscore sum → fractional bitscore per gene. This prevents large references from winning by sheer gene count. This would likely fix all 23 remaining failures without any database changes.

*Implementation: post-processing of the Kaptive TSV results file — no changes to Kaptive itself needed.*

#### Option B: Variable-region-only scoring

Remove conserved flanking/export genes (wza, wzb, wzc, galF, gnd, ugd) from the GenBank CDS features used by Kaptive, retaining only the variable biosynthetic region. Since these conserved genes are present at ~100% identity in all loci of the same KX type, they contribute equal bitscore to all references; removing them eliminates the shared background and leaves only the discriminating variable genes.

#### Expanded genome coverage

- Mine [AllTheBacteria](https://doi.org/10.1101/2024.03.08.584059) (~300–500K *E. coli* genomes) using [LexicMap](https://www.nature.com/articles/s41587-025-02812-8) or direct BLAST screening for *galF*/*gnd* flanking genes to discover novel G1/G4 K-locus types beyond the BSI isolate set

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
| **This database** | *E. coli* Group 1 & 4 (93 filtered loci, v0.3) | — |

## Citation

If you use this database, please cite:

- This repository
- Gladstone RA et al. (2024). *E. coli* capsule typing. [medRxiv preprint](https://www.medrxiv.org/content/10.1101/2024.11.22.24317484v1)
- The [EC-K-typing](https://github.com/rgladstone/EC-K-typing) Group 2 & 3 database

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
