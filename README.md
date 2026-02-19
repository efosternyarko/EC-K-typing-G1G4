# EC-K-typing Group 1 & 4

A reference database of *Escherichia coli* capsule (K-antigen) loci for **Group 1 and Group 4** capsule types.

This database complements the [EC-K-typing](https://github.com/rgladstone/EC-K-typing) Group 2 & 3 database by Rebecca Gladstone, enabling comprehensive capsule typing across all four *E. coli* capsule groups.

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

| Version | G1/G4 loci | All-groups loci | Source genomes | Notes |
|---------|-----------|-----------------|----------------|-------|
| v1.0 | 46 (KL300–KL343) | 136 | 222 (subset accessible) | Initial release |
| **v2.0** | **125 (KL300–KL423)** | **183** | **1,112 (all no-hit BSI isolates)** | **Current recommended version** |

### Reference loci (v2.0)

The v2.0 database contains **125 reference K-loci** (K24, K96, KL300–KL423) extracted from all 1,112 *E. coli* bloodstream infection isolates with no hit to the Group 2/3 database:

- **Size range (full set):** 7.0 – 59.1 kb
- **Total sequence (full set):** 4.14 Mb
- **Filtered set (≥ 30 kb):** 93 loci, 3.69 Mb — used in the combined all-groups database
- **Clustering threshold:** 95% nucleotide identity, 80% query coverage
- **KX types covered:** KX01, KX17, KX28, KX31, KX34, KX36, KX63, KX67, KX72, K24

### Files

#### v2.0 (current)

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v2.0.fasta` | All 125 reference locus sequences (FASTA) |
| `DB/EC-K-typing_group1and4_v2.0.gbk` | Annotated reference loci (GenBank, Kaptive-compatible) |
| `DB/EC-K-typing_group1and4_v2.0_filtered.fasta` | Filtered set of 93 loci ≥ 30 kb (FASTA) |
| `DB/EC-K-typing_all_groups_v2.0.gbk` | **Combined database** — all 4 capsule groups (183 loci, Kaptive-ready) |
| `DB/KL_G1G4_mapping.tsv` | KL nomenclature mapping for all 125 loci (KL name, KX origin, source assembly, length) |
| `DB/KL_G1G4_mapping_filtered.tsv` | Mapping for the 93-locus filtered set |
| `DB/cluster_info.tsv` | Clustering details (cluster members, representative sequences) |

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
5. **Representative selection:** Longest sequence per cluster

### Source data

The reference loci were derived from 6,673 *E. coli* bloodstream infection (BSI) isolates originally extracted from [EnteroBase](https://enterobase.warwick.ac.uk/):

1. **Group 2/3 screening:** All 6,673 BSI genomes were typed against the [Gladstone Group 2 & 3 database](https://github.com/rgladstone/EC-K-typing). **1,112 isolates** had no hit, indicating they likely carry Group 1 or Group 4 capsule loci.
2. **FastKaptive typing:** The 1,112 no-hit isolates were typed with [FastKaptive](https://github.com/rmostowy/fastKaptive), which assigned KX-type designations across 15 distinct types (KX36: 382, KX34: 329, KX17: 195, KX01: 73, KX67: 41, KX31: 29, KX72: 27, and others).
3. **Genome access:** All **1,112 genomes** were obtained for locus extraction in v2.0, resolving the EnteroBase accession mappings that limited v1.0 to 222 isolates.
4. **Locus extraction and clustering:** The extraction pipeline (see above) was applied to all 1,112 genomes, yielding **125 distinct K-locus clusters** (KL300–KL423, plus K24 and K96) covering 10 KX types. A filtered set of **93 loci ≥ 30 kb** is provided for use with Kaptive, excluding partial or fragmented extractions.

### Annotation

Gene prediction and annotation follows a two-stage approach:

1. **Gene prediction:** [pyrodigal](https://github.com/althonos/pyrodigal) in metagenomic mode identifies open reading frames across each locus
2. **Gene name transfer:** Predicted proteins are searched against the [Kaptive *Klebsiella* K-locus reference database](https://github.com/klebgenomics/Kaptive) via BLASTp (>=30% identity, >=50% coverage). Matching genes receive functional names from the *Klebsiella* reference (e.g., *wza*, *wzb*, *wzc*, *wzx*, *wzy*, glycosyltransferases). 73% of predicted CDS received gene name annotations.

GenBank records are formatted for [Kaptive](https://github.com/klebgenomics/Kaptive) compatibility, with each locus identified by a `source` feature containing `/note="K locus: KLxxx"`.

## Usage

### With Kaptive

A pre-built combined database (`DB/EC-K-typing_all_groups_v2.0.gbk`) covering all four capsule groups (183 loci) is included and ready for direct use with [Kaptive](https://github.com/klebgenomics/Kaptive):

```bash
kaptive assembly DB/EC-K-typing_all_groups_v2.0.gbk genome.fasta -o results.tsv
```

The combined database merges:
- **Group 2 & 3:** 90 loci (KL1–KL175) from [Gladstone et al.](https://github.com/rgladstone/EC-K-typing) — annotated with Bakta + Panaroo
- **Group 1 & 4:** 93 loci (K24, K96, KL300–KL423, filtered ≥ 30 kb) from this study — annotated with [pyrodigal](https://github.com/althonos/pyrodigal) (metagenomic mode), gene names transferred from the [Kaptive *Klebsiella* K-locus reference](https://github.com/klebgenomics/Kaptive) via BLASTp (73% of CDS annotated)

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

### v2.0 validation

The v2.0 database was validated by running Kaptive v3.1.0 on all 565 assemblies from which a locus was successfully extracted:

| Metric | Result |
|--------|--------|
| Self-typing (93 filtered loci) | **70/93 (75.3%)** |
| Self-typing (125 full set) | 70/125 (56.0%) |
| Assemblies typeable | **258/565 (45.7%)** |
| Loci utilised (filtered set) | 70/93 |

The 23 loci that do not self-type are predominantly KX01 loci that score slightly below KL302 (the longest KX01 representative at 41 kb). Reliable KL-level discrimination within a KX type requires systematic positional gene naming (planned for v3.0).

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

## Planned improvements

The following improvements are planned for future versions, drawing in part on approaches from [kTYPr](https://github.com/SushiLab/kTYPr) (Schwengers et al., 2025):

### Expanded genome coverage

- Mine [AllTheBacteria](https://doi.org/10.1101/2024.03.08.584059) (~300-500K *E. coli* genomes) using [LexicMap](https://www.nature.com/articles/s41587-025-02812-8) or direct BLAST screening for *galF*/*gnd* flanking genes to discover novel G1/G4 K-locus types

### Systematic gene naming

Adopt a positional gene naming strategy for locus-specific genes (as used by kTYPr):
- **Conserved genes:** functional names (e.g., *wza*, *wzb*, *wzc*, *wzx*, *wzy*, *galF*, *gnd*)
- **Locus-specific genes:** positional codes (e.g., *KL300_7*, *KL300_8*)
- **Shared variable genes:** compound names reflecting shared usage across types (e.g., *KL300_KL305_7*)
- Protein families defined by clustering at 90% identity/coverage, with one HMM profile per family for robust detection

### Enhanced scoring and reporting

- **Multi-metric scoring** per K-type assignment (as in kTYPr): number of reference genes found, accumulated bitscore, fractional bitscore, and locus completeness
- **Conserved locus integrity check:** separate scoring for conserved flanking/export genes (*wza*, *wzc*, *galF*, *gnd*) vs. variable biosynthetic genes
- **Locus comparison visualisation:** [clinker](https://github.com/gamcil/clinker)-based HTML visualisation for pairwise comparison of K-locus architectures

## Related work

| Tool/Database | Scope | Reference |
|---------------|-------|-----------|
| [EC-K-typing](https://github.com/rgladstone/EC-K-typing) | *E. coli* Group 2 & 3 (90 loci) | [Gladstone et al. 2024](https://www.medrxiv.org/content/10.1101/2024.11.22.24317484v1) |
| [kTYPr](https://github.com/SushiLab/kTYPr) | *E. coli* Group 2 & 3 (85 loci, HMM-based) | [Schwengers et al. 2025](https://www.biorxiv.org/content/10.1101/2025.08.07.669119v1) |
| [Kaptive](https://github.com/klebgenomics/Kaptive) | *Klebsiella* and *E. coli* K/O typing | [Lam et al. 2022](https://doi.org/10.1099/mgen.0.000800) |
| [FastKaptive](https://github.com/rmostowy/fastKaptive) | Fast K-locus pre-screening | Mostowy et al. |
| **This database** | *E. coli* Group 1 & 4 (125 loci, v2.0) | — |

## Citation

If you use this database, please cite:

- This repository
- Gladstone RA et al. (2024). *E. coli* capsule typing. [medRxiv preprint](https://www.medrxiv.org/content/10.1101/2024.11.22.24317484v1)
- The [EC-K-typing](https://github.com/rgladstone/EC-K-typing) Group 2 & 3 database

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
