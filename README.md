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

### Reference loci

The database contains **46 reference K-loci** (K24, K96, KL300-KL343) extracted from bloodstream infection *E. coli* isolates:

- **Size range:** 33.4 - 46.1 kb
- **Total sequence:** 1.77 Mb
- **Clustering threshold:** 95% nucleotide identity, 80% query coverage
- **KX types covered:** 19 of 20 FastKaptive-assigned types

### Files

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v1.0.fasta` | Reference locus sequences (FASTA) |
| `DB/EC-K-typing_group1and4_v1.0.gbk` | Annotated reference loci (GenBank, Kaptive-compatible) |
| `DB/KL_G1G4_mapping.tsv` | KL nomenclature mapping (KL name, KX origin, source assembly, length) |
| `DB/cluster_info.tsv` | Clustering details (cluster members, representative sequences) |
| `scripts/build_G1G4_db.py` | Pipeline script for locus extraction and database construction |
| `scripts/annotate_loci.py` | Annotation script (pyrodigal gene prediction + BLAST annotation) |
| `flanking_genes/flanking_genes.fasta` | Flanking gene marker sequences used for locus detection |

### Nomenclature

| KL designation | Description |
|----------------|-------------|
| K24, K96 | Direct K-type calls (known capsule types) |
| KL300-KL343 | Novel locus types identified in this study |

## Methods

### Locus extraction pipeline

1. **Flanking gene detection:** BLAST search for *galF*, *gnd*, *ugd*, *wza*, and *wzc* from *E. coli* K-12 MG1655 against genome assemblies
2. **Locus boundary definition:** The *cps* region between *galF* and *gnd* (including *wza-wzb-wzc* upstream) with 500 bp flanking sequence
3. **Fragmented assembly handling:** For genomes where *galF* and *gnd* are on different contigs, sequences are extracted from both contigs and concatenated with an N-spacer
4. **Clustering:** All-vs-all BLAST at 95% identity and 80% query coverage, greedy clustering by sequence length
5. **Representative selection:** Longest sequence per cluster

### Source data

- 219 genome assemblies from the *E. coli* bloodstream infection collection (207 downloaded from ENA + 12 from EnteroBase)
- Untypeable isolates from [Gladstone et al.](https://www.medrxiv.org/content/10.1101/2024.11.22.24317484v1) that did not match the Group 2 & 3 database
- Initial K-locus type assignments from [FastKaptive](https://github.com/rmostowy/fastKaptive)

## Usage

### With Kaptive

The GenBank file (`DB/EC-K-typing_group1and4_v1.0.gbk`) is ready for use with [Kaptive](https://github.com/klebgenomics/Kaptive). To create a combined database covering all four capsule groups:

```bash
cat EC-K-typing_group2and3_v3.0.0.gbk DB/EC-K-typing_group1and4_v1.0.gbk \
  > EC-K-typing_all_groups_v1.0.gbk

kaptive assembly -k EC-K-typing_all_groups_v1.0.gbk -a genome.fasta
```

Gene predictions were generated using [pyrodigal](https://github.com/althonos/pyrodigal) (metagenomic mode) with known capsule pathway genes (*galF*, *gnd*, *ugd*, *wza*, *wzc*) annotated via tBLASTn against *E. coli* K-12 MG1655 reference sequences.

### With BLAST

For direct BLAST-based typing:

```bash
makeblastdb -in DB/EC-K-typing_group1and4_v1.0.fasta -dbtype nucl

blastn -query genome.fasta \
  -db DB/EC-K-typing_group1and4_v1.0.fasta \
  -outfmt '6 qseqid sseqid pident qcovs length' \
  -max_target_seqs 5
```

## Citation

If you use this database, please cite:

- This repository
- Gladstone RA et al. (2024). *E. coli* capsule typing. [medRxiv preprint](https://www.medrxiv.org/content/10.1101/2024.11.22.24317484v1)
- The [EC-K-typing](https://github.com/rgladstone/EC-K-typing) Group 2 & 3 database

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
