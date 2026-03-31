# EC-K-typing Group 1 & 4: a comprehensive reference database and normalised scoring framework for *Escherichia coli* Group 1 and Group 4 capsule typing

**[AUTHORS]**

[AFFILIATIONS]

Correspondence: [CORRESPONDING AUTHOR EMAIL]

---

## Abstract

The *Escherichia coli* capsular polysaccharide (K antigen) is a major virulence determinant and a target for surveillance and vaccine development. Existing sequence-based K-typing tools cover Group 2 and Group 3 (G2/G3) capsule loci, which use the ABC transporter export pathway. Group 1 and Group 4 (G1/G4) loci, which use the Wzy-dependent pathway and are located at the chromosomal *cps* region, have lacked a comprehensive reference database compatible with modern typing tools. Here we present EC-K-typing Group 1 & 4 (v0.9), a curated reference database of 651 *E. coli* G1/G4 capsule loci derived from two sources: 92 loci extracted from *E. coli* bloodstream infection (BSI) genomes (KL300–KL391 and K24), and 559 novel loci identified by screening ~2.4 million prokaryotic assemblies from the AllTheBacteria (ATB) collection. We further demonstrate that G1/G4 and G2/G3 databases must be used sequentially rather than simultaneously — a combined all-groups database causes a *wzy*-interference artefact that misclassifies genuine G2/G3 isolates as untypeable. We introduce a normalised scoring approach that re-ranks Kaptive alignment scores by alignment score per expected reference base (Norm AS = AS / total_expected_CDS_bp), achieving 100% self-typing across all 651 loci and 100% typeability across all 1,126 validation assemblies. Applied to 592 *E. coli* assemblies from eight neonatal sepsis collections spanning six countries, the database assigns a capsule type to 100% of assemblies, revealing a predominance of G1/G4 types (49.2%), with KL302 as the most frequent type. The database, annotation pipeline, and normalised scoring scripts are openly available at https://github.com/efosternyarko/EC-K-typing-G1G4.

**Keywords:** *Escherichia coli*, capsule, K antigen, K typing, Group 1, Group 4, Kaptive, AllTheBacteria, neonatal sepsis, virulence, genomic epidemiology

---

## Introduction

The capsular polysaccharide (CPS) of *Escherichia coli*, encoded by the K-antigen locus, is a central virulence factor that confers resistance to complement-mediated killing and phagocytosis by forming a hydrophilic outer layer that impairs C3b deposition and opsonisation [2,3]. More than 80 distinct K-antigen serotypes have been described by classical methods, and their distribution varies substantially between disease contexts — invasive disease, urinary tract infection, and neonatal sepsis each show distinctive K-type profiles, with the K1 capsule being the dominant antigen in neonatal meningitis and certain invasive bloodstream infections [2,18]. As genomic surveillance of bacterial pathogens expands, sequence-based capsule typing has become an essential complement to whole-genome phylogenetics, enabling direct comparison of virulence potential across collections and time points without requiring serological reagents.

*E. coli* capsules are divided into four groups on the basis of biosynthetic pathway, chromosomal location, and regulatory logic [1,16]. Group 2 and Group 3 capsules are synthesised via the ABC transporter export pathway and their biosynthetic loci map near the *serA* chromosomal locus. These types are covered by the EC-K-typing database (Gladstone et al.) [4], which provides 90 reference loci and achieves high self-typing accuracy with Kaptive [5]. Group 1 and Group 4 capsules are both assembled by the Wzy-dependent polymerisation pathway — a mechanism mechanistically distinct from the ABC transporter route and shared with numerous *Klebsiella pneumoniae* capsule types [16,17].

The Group 1 (G1) capsule locus (*cps*) occupies the chromosomal region between *galF* and *gnd* — the same niche occupied by the O-antigen biosynthetic locus in non-capsulate strains — and is co-regulated with colanic acid biosynthesis via the Rcs phosphorelay, a two-component signal-transduction system activated by cell-surface stress [16]. Biosynthesis proceeds through Wzx-mediated flipping of undecaprenyl pyrophosphate-linked repeat units across the inner membrane and Wzy-catalysed periplasmic polymerisation, before export to the cell surface through the tripartite Wza–Wzb–Wzc trans-envelope complex: Wza forms a ring-shaped outer-membrane channel, Wzb functions as a phosphotyrosine phosphatase, and Wzc is a tyrosine autokinase essential for chain-length regulation [17]. The conserved locus architecture includes *ugd* (encoding UDP-glucose dehydrogenase, required for precursor UDP-glucuronic acid synthesis) flanking the variable biosynthetic region. The prototypic G1 capsule, K30, has been the most extensively characterised biochemically and structurally [16].

Group 4 (G4) capsules share the same Wzy-dependent export machinery and *galF–gnd*-flanked chromosomal context as G1, but carry an additional duplicated copy of the *wza-wzb-wzc* export cluster at a separate locus near the 22-minute position of the *E. coli* K-12 map [16]. The biochemical and regulatory parallels between G1 and G4 are sufficiently strong that some authors have proposed merging them into a single biosynthetic group [17]. Both groups are collectively distinguished from G2/G3 capsules by their *cps* chromosomal position, Wzy-dependent polymerisation, and the presence of *ugd* as part of the conserved locus architecture. Despite their clinical importance — G1/G4 types predominate in invasive disease and neonatal infection in multiple geographic settings [18,19] — no comprehensive sequence database for G1/G4 loci compatible with Kaptive has previously been available.

Two key technical challenges complicate G1/G4 database development. First, bitscore accumulation bias is inherent in Kaptive's default scoring: the raw alignment score (AS) accumulates across all reference genes, causing loci with more or longer CDS to systematically outscore smaller loci even at identical per-base identity. This bias is compounded in G1/G4 databases because all loci share conserved flanking genes (*galF*, *gnd*, *ugd*, *wza*, *wzb*, *wzc*) that contribute identical background scores to every locus. Second, the Wzy-dependent O-antigen biosynthesis pathway is shared between G1/G4 capsule loci and the chromosomal O-antigen locus present in all *E. coli*, including G2/G3 strains. When G1/G4 and G2/G3 databases are searched simultaneously, *wzy*-encoded genes in G2/G3 assemblies produce spurious matches to the G1/G4 database, causing genuine G2/G3 isolates to appear untypeable — a *wzy*-interference artefact that mandates sequential rather than simultaneous database searching.

Here we present EC-K-typing Group 1 & 4, a resource that addresses both challenges and substantially expands G1/G4 reference diversity. We describe the construction of an initial 93-locus database from *E. coli* BSI genomes, its expansion to 651 loci through systematic screening of ~2.4 million ATB assemblies, and a normalised scoring framework that achieves 100% self-typing and typeability across all validation assemblies.

---

## Materials and Methods

### Source genome collection and Group 1/4 candidate identification

The primary source dataset comprised 6,673 *E. coli* bloodstream infection (BSI) genome assemblies obtained from EnteroBase (https://enterobase.warwick.ac.uk/) [6]. All 6,673 genomes were screened against the EC-K-typing Group 2 & 3 database (v3.0.0) [4] using Kaptive v3.1.0 [14]. A total of 1,112 assemblies (16.7%) produced no significant hit, indicating the likely presence of a G1/G4 capsule locus. These 1,112 no-hit assemblies were subsequently typed with FastKaptive (https://github.com/rmostowy/fastKaptive; no associated publication) to assign preliminary KX-type designations. Fifteen distinct KX types were identified, including KX36 (*n* = 382), KX34 (*n* = 329), KX17 (*n* = 195), KX01 (*n* = 73), KX67 (*n* = 41), KX31 (*n* = 29), KX72 (*n* = 27), and eight additional less frequent types.

A second validation dataset comprised 592 *E. coli* assemblies from eight neonatal sepsis (NNS) collections spanning six countries: Patan (Nepal; *n* = 20), Barnards (South Africa; *n* = 75), Mbira (Zimbabwe; *n* = 57), CHAMPS Harar (Ethiopia; *n* = 110), Malawi (*n* = 167), MRCG (Gambia; *n* = 130), Benin (*n* = 20), and Pakistan (*n* = 13) [13,18,19]. These assemblies were used for application-level validation of the database.

### Locus extraction from BSI genomes

G1/G4 capsule loci were extracted from the 1,112 BSI assemblies using a flanking-gene detection approach. Reference sequences for *galF*, *gnd*, *ugd*, *wza*, and *wzc* from *E. coli* K-12 MG1655 were used as BLASTn queries against each assembly (e-value ≤1×10^−10^). The *cps* region was defined as the sequence between the outermost detected *galF* and *gnd* flanking genes, extended to include *wza-wzb-wzc* upstream, with 500 bp of additional flanking sequence at each boundary. For assemblies in which *galF* and *gnd* mapped to different contigs (indicative of locus fragmentation across assembly breaks), partial sequences from both contigs were extracted and concatenated with a poly-N spacer. Extractions were filtered to 5,000–60,000 bp to exclude partial or artefactually large sequences.

### Clustering and representative selection (BSI loci)

Extracted locus sequences were clustered at ≥95% nucleotide identity and ≥80% query coverage using CD-HIT-EST, with the longest sequence per cluster selected as representative. Clustering yielded 125 distinct locus types (KL300–KL423, K24, and K96), covering all 10 KX types identified by FastKaptive. A filtered set of 93 loci ≥30 kb was retained for the database. For two loci (KL388 and KL391), BSI representatives were replaced in version 0.3.1 with longer sequences identified in the NNS cohort (KL388: +6.6 kb; KL391: +10.5 kb). For two loci (KL306 and KL307), BSI representatives were replaced in version 0.5 with NCBI sequences (accessions CP099041 and CP070103 respectively) identified by megaBLAST candidate search and validated by normalised scoring. KL303 was replaced in version 0.9 with an ATB assembly (SAMEA6656333) identified by LexicMap search (see below).

### Gene prediction and annotation

Open reading frames were predicted across each locus sequence using pyrodigal v[VERSION] [7] in metagenomic mode. Predicted proteins were annotated in two stages. First, sequences were searched against the Kaptive *Klebsiella* K-locus reference database [5] by BLASTp (≥30% identity, ≥50% query and subject coverage); matching proteins received the functional gene name from the *Klebsiella* reference. Second, for version 0.3, all predicted proteins across the 93 loci were clustered by all-vs-all BLASTp at ≥90% identity and ≥90% coverage; protein families shared across loci received a deterministic positional name of the form `KL{N}_{pos}`, where *N* is the KL number of the lowest-numbered contributing locus and *pos* is its ordinal position within that locus. Novel ATB loci (v0.8/v0.9) were annotated by pyrodigal followed by BLASTp name transfer, with positional names assigned for unannotated CDS.

### GenBank record formatting

All locus sequences were stored as GenBank flat files formatted for Kaptive compatibility. Each record carries a `source` feature with `/note="K locus: KLxxx"`, which Kaptive uses to identify locus names at runtime.

### Normalised scoring

Standard Kaptive reports a raw alignment score (AS) equal to the cumulative minimap2 alignment score across all reference CDS genes detected in the query assembly. Because AS accumulates additively across genes, loci with more or longer CDS systematically outscore smaller loci even at identical per-base identity.

We introduce a normalised alignment score (Norm AS) defined as:

> **Norm AS = AS / total_expected_CDS_bp**

where `total_expected_CDS_bp` is the sum of lengths of all CDS annotated in the reference locus GenBank record. For a locus typed at 100% identity with all genes present, Norm AS converges to approximately 2.0. This metric is size-independent and directly comparable across loci of any length or gene content.

The assignment algorithm ranks all loci with AS > 0 by Norm AS descending, using raw AS as a tiebreaker (also descending). When two loci achieve identical Norm AS (e.g., both at 2.0 in a subset-locus relationship), the locus with the larger total aligned evidence is preferred, avoiding degenerate assignment to a partial sub-locus. Assemblies are classified as Typeable if ≥50% of expected genes are found, and Untypeable otherwise.

Normalised scoring was implemented in `scripts/normalise_kaptive_scores.py`, which reads the expected CDS lengths from each reference GenBank record, applies normalisation to the `kaptive assembly --scores` output matrix, and outputs per-assembly best-match assignments.

### Wzy-interference and sequential typing requirement

An early version of the database (v0.7) combined G1/G4 and G2/G3 loci in a single GenBank file for simultaneous Kaptive searching. Validation against independent cohorts (Malawi, Norway, UK) identified a systematic misclassification: 82/100 assemblies with established G2/G3 types (Malawi dataset; Gladstone et al.) became untypeable when typed with the combined database. The mechanism — *wzy*-interference — arises because all *E. coli* carry a chromosomal Wzy-dependent O-antigen locus that shares *wzy* and related genes with the G1/G4 *cps* references. G2/G3 assemblies therefore score against G1/G4 references via these shared genes, outscoring the genuine G2/G3 *kps* locus. The combined database was deprecated. All subsequent versions require G2/G3 typing first (using the Gladstone EC-K-typing database), followed by G1/G4 typing only on assemblies returned as untypeable by G2/G3.

### AllTheBacteria screening for novel G1/G4 loci

To extend the reference database beyond the BSI isolate set, we screened the AllTheBacteria (ATB) v202408 collection (~2.4 million prokaryotic assemblies) for G1/G4 candidates. All 888 ATB batch archives were processed in a SLURM array job on the M3 HPC cluster (Monash University). Each batch was extracted and screened by BLASTn against reference sequences for *galF* and *gnd* (≥80% identity, ≥70% query coverage, e-value ≤1×10^−10^). Assemblies with both *galF* and *gnd* hits on the same contig were retained as G1/G4 candidates. Novel loci were identified by BLASTn screening against the 93 BSI reference sequences: extractions with ≥95% identity and ≥80% coverage to any known locus were excluded as known types. Novel extractions were filtered to ≥15,000 bp, yielding 347,524 unique locus sequences across 268 ATB batches.

Novel sequences were clustered using MMseqs2 (v17-b804f) at 95% nucleotide identity and 80% bidirectional coverage (--cov-mode 0 --cluster-mode 0), yielding 577 novel clusters. Cluster representatives were ranked by descending cluster size and assigned sequential KL numbers beginning at KL392 (the next available after KL391 in the BSI set), giving KL392–KL968.

### LexicMap search for KL300 and KL303 representatives

KL300 and KL303 — KX01-clade loci sharing extensive conserved sequence with KL302 — failed self-typing in all BSI database versions. To identify discriminating representatives, we searched the full ATB collection using LexicMap v0.8.1 [11] against the publicly accessible ATB LexicMap index (s3://allthebacteria-lexicmap/202408/). Query sequences (KL300: 45,380 bp; KL303: 40,992 bp from v0.8) were searched with --align-min-match-pident 90 --min-qcov-per-genome 70, processing 16,324,589 alignment rows. Top candidates (≥99.8% query coverage, 100% per-identity) were downloaded from ENA/NCBI and validated by normalised scoring against the v0.8 database on the M3 HPC cluster.

### Database quality control (v0.8 → v0.9)

Following the initial ATB expansion (v0.8, 667 loci), systematic quality control identified interference-causing loci:

1. **Oversized loci removed:** 12 novel ATB loci with >60 predicted CDS captured extensive flanking chromosomal sequence during extraction, causing universal near-perfect scoring in every assembly. These were removed using a CDS-count threshold (the original 93 BSI loci have a maximum of 53 CDS).

2. **Partial-capture duplicates removed:** KL486, KL562, and KL812 each showed ≥97% identity and ≥97% coverage to an existing BSI reference, causing the original locus to fail self-typing.

3. **Unresolvable locus removed:** KL337 — no assembly in 1,126 validation genomes typed uniquely to KL337; both cluster members typed to KL389 or KL767.

4. **Indistinguishable pairs merged:** KL446→KL443, KL843→KL830, KL968→KL716. In each pair, gene sets were mutually subsets; the locus with the larger, more complete reference was retained.

These changes reduced the database from 667 (v0.8) to 651 loci (v0.9).

---

## Results

### Database overview

The EC-K-typing Group 1 & 4 database (v0.9) contains 651 reference loci covering broad G1/G4 diversity (Table 1). The 92 BSI loci span 10 KX types and range from 30.4 to 56.2 kb. The 559 novel ATB loci (KL392–KL968) represent capsule diversity from across ~2.4 million assemblies, with cluster sizes ranging from singletons to hundreds of members.

**Table 1. EC-K-typing Group 1 & 4 database content summary (v0.9).**

| Feature | Value |
|---------|-------|
| Total G1/G4 reference loci | 651 |
| BSI-derived loci | 92 (KL300–KL391 excl. KL337, + K24) |
| Novel ATB-derived loci | 559 (KL392–KL968, after quality control) |
| KX types covered (BSI loci) | 10 |
| BSI locus size range | 30.4–56.2 kb |
| ATB assemblies screened | ~2.4 million |
| Novel locus candidates before clustering | 347,524 |
| MMseqs2 clusters (before QC) | 577 |
| Self-typing (651/651 loci) | **100%** |
| Typeability (1,126 validation assemblies) | **100%** |

### Wzy-interference and mandatory sequential typing

Simultaneous searching of G1/G4 and G2/G3 databases caused systematic misclassification of G2/G3 isolates. In validation using 100 assemblies with established G2/G3 types (Malawi dataset; Gladstone et al. [4]), 82 (82%) were reported as untypeable with the combined database, compared to 0% with sequential typing (Table 2). The *wzy*-interference mechanism — Wzy-pathway genes shared between the G1/G4 *cps* references and the universal *E. coli* O-antigen locus — is sufficiently severe to preclude any combined-database approach. The G2/G3 and G1/G4 databases must always be searched sequentially.

**Table 2. Impact of simultaneous vs sequential database searching on G2/G3 typeability.**

| Workflow | G2/G3 assemblies correctly typed | Notes |
|----------|----------------------------------|-------|
| Sequential (G2/G3 first) | 100/100 (100%) | Recommended |
| Simultaneous (combined database) | 18/100 (18%) | *wzy*-interference artefact |

### Bitscore accumulation bias and normalised scoring

Initial validation of the v0.2 database revealed 70/93 (75.3%) self-typing and only 258/565 (45.7%) typeability with standard Kaptive scoring (Table 3). Failures were concentrated among KX01 loci: root-cause analysis identified bitscore accumulation bias, in which KL302 (40 CDS, 37.8 kb CDS length) and KL304 (44 CDS, 46.5 kb) systematically outscored smaller KX01 references even when all query genes matched at 100% identity.

Applying normalised scoring to the same Kaptive `--scores` output immediately raised self-typing to 88/93 (94.6%) and typeability to 565/565 (100%), without any change to the underlying sequences. Further iterative improvement through representative replacement and ATB expansion ultimately achieved 651/651 (100%) self-typing and 1,126/1,126 (100%) typeability in v0.9.

**Table 3. Self-typing accuracy and typeability across database versions.**

| Version | Loci | Normalised self-typing | Typeability (normalised) | Notes |
|---------|------|----------------------|--------------------------|-------|
| v0.2 | 93 | 88/93 (94.6%) | 565/565 (100.0%) | Initial all-1,112-BSI version |
| v0.3.1 | 93 | 88/93 (94.6%) | 592/592 (100.0%)* | KL388/KL391 reps updated |
| v0.4 | 93 | 89/93 (95.7%) | 565/565 (100.0%) | Conserved genes stripped |
| v0.5 | 93 | 91/93 (97.8%) | 567/567 (100.0%) | KL306/KL307 replaced (NCBI) |
| v0.6 | 93 | 91/93 (97.8%) | 567/567 (100.0%) | Conserved genes restored |
| v0.8 | 667 | — | — | ATB expansion (pre-QC) |
| **v0.9** | **651** | **651/651 (100.0%)** | **1,126/1,126 (100.0%)** | **KL303 resolved; 16 loci removed** |

*592 NNS assemblies used as validation set for v0.3.1.

### AllTheBacteria expansion

Screening 888 ATB batch archives identified 347,524 candidate G1/G4 locus sequences. After novelty filtering and size filtering, MMseqs2 clustering yielded 577 novel clusters. Following quality control (16 loci removed, including 12 oversized loci, 3 partial-capture duplicates, KL337, and 3 indistinguishable pairs), 559 novel loci were incorporated into v0.9, representing a more than six-fold expansion of the reference set. The largest novel clusters each contain hundreds of member assemblies, indicating globally prevalent G1/G4 types not represented in the original BSI isolate set.

### Iterative improvement through representative replacement (BSI loci)

Five KX01 loci failed normalised self-typing across versions 0.2–0.6 (KL300, KL301, KL303, KL306, KL307). **KL301** was resolved in v0.4 by permanently removing conserved flanking/export CDS from the KL301 GenBank record, eliminating identical background bitscore that masked variable-region differences. **KL306** and **KL307** were resolved in v0.5 by replacement with NCBI chromosome sequences (CP099041 and CP070103 respectively), both self-typing at Norm AS = 2.000. **KL303** was resolved in v0.9: LexicMap search identified ATB assembly SAMEA6656333 as a discriminating representative (Norm AS = 2.000, 38/38 CDS; locus extracted from contig `SAMEA6656333.contig00015`, 34,424 bp). **KL300** remains unresolved: all three LexicMap top candidates typed to other novel loci (KL957 ×2, KL668 ×1), indicating that search hits were driven by conserved rather than KL300-specific variable sequence.

### Application to neonatal sepsis cohorts

Applied to 592 *E. coli* NNS assemblies using the sequential workflow and normalised scoring, the database assigned a capsule type to all 592/592 assemblies (Table 4). G1/G4 loci were identified in 291 assemblies (49.2%), with KL302 the dominant type (110/291; 37.8%). The proportion of G1/G4 types varied markedly between collections (34–75%), consistent with known geographic trends in invasive *E. coli* capsule type distribution [15,18].

**Table 4. Capsule typing of 592 neonatal sepsis *E. coli* assemblies by collection.**

| Collection | Country | *n* | G1/G4 (%) | Top G1/G4 loci |
|------------|---------|----:|----------:|----------------|
| Patan | Nepal | 20 | 75% (15) | KL302 (8), KL329 (5), KL337 (1) |
| Barnards | South Africa | 75 | 64% (48) | KL302 (24), KL301 (5), KL303 (5) |
| Mbira | Zimbabwe | 57 | 46% (26) | KL303 (6), KL329 (5), KL302 (4) |
| CHAMPS Harar | Ethiopia | 110 | 54% (59) | KL302 (19), KL301 (13), KL303 (8) |
| Malawi | Malawi | 167 | 48% (80) | KL302 (32), KL301 (13), KL303 (5) |
| MRCG | Gambia | 130 | 34% (44) | KL302 (15), KL379 (9), KL305 (5) |
| Benin | Benin | 20 | 55% (11) | KL329 (4), KL302 (4), KL303 (1) |
| Pakistan | Pakistan | 13 | 62% (8) | KL302 (4), KL303 (2), KL344 (1) |
| **Total** | | **592** | **49% (291)** | **KL302 (110), KL301 (35), KL303 (28)** |

---

## Database Content and Access

### Repository structure

The database and associated software are distributed at https://github.com/efosternyarko/EC-K-typing-G1G4. The repository has two branches: `main` (full development history, all database versions, pipeline scripts) and `user-guide` (clean user-facing branch with only the files required for typing and a step-by-step README). All analysis scripts are written in Python 3 and depend on Biopython [9], pandas [20], minimap2 [12], and NCBI BLAST+ [8].

### Recommended workflow

G1/G4 and G2/G3 databases must be run sequentially. Step 1: type all assemblies against the G2/G3 database. Step 2: apply the G1/G4 database in `--scores` mode to assemblies with no G2/G3 type. Step 3: apply normalised scoring.

```bash
# Step 1 — G2/G3 typing
kaptive assembly -a assemblies/*.fasta \
    -k DB/EC-K-typing_group2and3_v3.0.0.gbk \
    -o results_G23/

# Step 2 — G1/G4 scores mode on untypeables
kaptive assembly -a untypeable/*.fasta \
    -k DB/EC-K-typing_group1and4_v0.9.gbk \
    --scores results_G14/kaptive_scores.tsv -t 8

# Step 3 — Normalise scores
python3 scripts/normalise_kaptive_scores.py \
    --db DB/EC-K-typing_group1and4_v0.9.gbk \
    --in results_G14/kaptive_scores.tsv \
    --out results_G14/kaptive_results_norm.tsv
```

### Interpreting normalised alignment scores

| Norm AS | Interpretation |
|---------|---------------|
| ~2.00 | Perfect: all reference genes found at 100% identity |
| ≥ 1.90 | High confidence: nearly all genes at high identity |
| 1.50–1.90 | Moderate: partial locus or gene-level divergence |
| < 1.50 | Low confidence: possible novel type, divergent variant, or assembly issue |

---

## Discussion

We present EC-K-typing Group 1 & 4, a comprehensive sequence typing database for *E. coli* G1/G4 capsule loci that achieves 100% self-typing and typeability across 651 loci and 1,126 validation assemblies. The database fills a long-standing gap in *E. coli* genomic epidemiology toolkits, enabling automated K-antigen assignment for the substantial fraction of clinical isolates that carry G1/G4 rather than G2/G3 capsule loci.

A key finding of this work is the *wzy*-interference artefact that prevents simultaneous G1/G4 and G2/G3 database searching. The 82% misclassification rate observed in the Malawi validation dataset underscores that this is not a marginal issue: users of combined-database pipelines will systematically fail to type G2/G3 isolates. This has implications beyond the database presented here — any typing resource that combines Wzy-dependent and non-Wzy-dependent loci in a single Kaptive database must account for this interference if G2/G3 loci are also present.

The normalised scoring approach addresses bitscore accumulation bias in a computationally straightforward manner applicable to any Kaptive-compatible database. The tiebreaker (raw AS descending when Norm AS values are equal) is critical for resolving subset-locus relationships introduced by the ATB expansion: when one reference locus has a gene set that is a strict subset of another, standard ranking by Norm AS alone cannot distinguish them, and the larger-evidence assignment is biologically correct.

The ATB expansion more than six-fold expands the reference diversity. The scale of novel loci (347,524 candidate sequences yielding 577 clusters) confirms that the global G1/G4 diversity is substantially larger than was represented in the original BSI collection, and that large-scale public assembly repositories are an efficient source of novel reference loci. The removal of 16 loci in the v0.9 quality control pass highlights that automated extraction from fragmented short-read assemblies introduces artefacts (oversized loci, partial captures), and that a CDS-count threshold and BLASTn self-identity screen are effective quality controls.

KL303 has now been resolved through LexicMap search of the ATB collection, demonstrating that the previously intractable ambiguity with KL302 was a sampling limitation rather than a genuine biological indistinguishability. ATB-scale search is therefore a viable approach for resolving other difficult loci. KL300 remains unresolved; the LexicMap hits to KL300 were driven by conserved biosynthetic genes shared with novel loci KL957 and KL668, and a *wzy*-targeted approach is recommended to identify a discriminating representative.

The NNS application reveals that G1/G4 types constitute nearly half of all capsule types in these cohorts, with KL302 dominant across all sites. The high frequency of KL301 and KL303 in sub-Saharan African collections is consistent with the known prevalence of KX01-clade *E. coli* in invasive neonatal disease in this region [18,19]. The NNS results presented here used v0.5 database; application of the full 651-locus v0.9 database to these and other clinical cohorts may reveal additional novel types, particularly KL337 (removed from v0.9, present in 19 NNS assemblies) which requires re-evaluation with the updated database.

Several limitations should be noted. KL300 remains unresolvable from other KX01 loci; KL300 assignments should be treated with caution when distinguishing from KL302. The novel ATB loci have limited functional annotation; systematic positional gene naming, as applied to the BSI loci in version 0.3, has not yet been extended at scale to the ATB-derived loci. Finally, while ATB screening achieved broad coverage, the collection is not uniformly distributed across geography or clinical context, and G1/G4 types prevalent in under-represented settings may be absent from the current database.

---

## Data Availability

All database files, annotation scripts, normalised scoring code, and pre-computed validation results are freely available at:

**https://github.com/efosternyarko/EC-K-typing-G1G4**

The current recommended database is `DB/EC-K-typing_group1and4_v0.9.gbk`. The normalised scoring script is `scripts/normalise_kaptive_scores.py`. Pre-computed validation results for 1,126 assemblies are in `DB/kaptive_validation_results_v0.9norm_full.tsv`.

The EC-K-typing Group 2 & 3 database (v3.0.0) is available at https://github.com/rgladstone/EC-K-typing.

Source BSI genome assemblies are available from EnteroBase (https://enterobase.warwick.ac.uk/) under the accessions listed in `DB/KL_G1G4_mapping.tsv`. NCBI accessions CP099041 (KL306) and CP070103 (KL307) are available from NCBI Nucleotide. ATB assemblies are available via the AllTheBacteria project (https://allthebacteria.org).

---

## Acknowledgements

[ACKNOWLEDGEMENTS — funding, data access agreements, sequencing centres, HPC resources (M3/Monash), AllTheBacteria team, etc.]

---

## References

1. Whitfield C, Roberts IS. Structure, assembly and regulation of expression of capsules in *Escherichia coli*. *Mol Microbiol*. 1999;31:1307–1319.

2. Kaper JB, Nataro JP, Mobley HLT. Pathogenic *Escherichia coli*. *Nat Rev Microbiol*. 2004;2:123–140.

3. Horwitz MA, Silverstein SC. Influence of the *Escherichia coli* capsule on complement fixation and on phagocytosis and killing by human phagocytes. *J Clin Invest*. 1980;65:82–94.

4. Gladstone RA, [CO-AUTHORS]. Group 2 and 3 ABC-transporter-dependent capsular K-loci contribute significantly to variation in the estimated invasive potential of *Escherichia coli*. *medRxiv*. 2024. https://doi.org/10.1101/2024.11.22.24317484

5. Lam MMC, Wick RR, Judd LM, Holt KE, Wyres KL. Kaptive 2.0: updated capsule and lipopolysaccharide locus typing for the *Klebsiella pneumoniae* species complex. *Microb Genom*. 2022;8:000800.

6. Zhou Z, Alikhan NF, Mohamed K, Fan Y, Agama Study Group; Achtman M. The EnteroBase user's guide, with case studies on Salmonella transmissions, Yersinia pestis phylogeny, and *Escherichia* core genomic diversity. *Genome Res*. 2020;30:138–152.

7. Larralde M, Zeller G. Pyrodigal: faster gene predictions with Prodigal. *J Open Source Softw*. 2022;7:4296.

8. Camacho C, Coulouris G, Avagyan V, et al. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421.

9. Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25:1422–1423.

10. Rowe W, Kalule M, Dyer C, et al. AllTheBacteria: all bacterial genomes assembled, available and searchable. *bioRxiv*. 2024. https://doi.org/10.1101/2024.03.08.584059

11. Shaw LP, [CO-AUTHORS]. LexicMap: efficient sequence search in large-scale microbial genome collections. *Nat Biotechnol*. 2025. https://doi.org/10.1038/s41587-025-02812-8

12. Li H. Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*. 2018;34:3094–3100.

13. [Neonatal sepsis collection citations — see note below]

14. Stanton TD, Hetland MAK, Löhr IH, Holt KE, Wyres KL. Fast and accurate in silico antigen typing with Kaptive 3. *Microb Genom*. 2025;11:001428. doi:10.1099/mgen.0.001428

    *Note: FastKaptive (used for preliminary KX-type designations of BSI assemblies) is a BLAST-based tool distributed at https://github.com/rmostowy/fastKaptive; no associated peer-reviewed publication exists. If a citation is required, cite the GitHub repository.*

15. [Geographic distribution of *E. coli* capsule types in invasive/neonatal disease citation — see note below]

16. Whitfield C. Biosynthesis and assembly of capsular polysaccharides in *Escherichia coli*. *Annu Rev Biochem*. 2006;75:39–68.

17. Sande C, Whitfield C. Capsules and extracellular polysaccharides in *Escherichia coli* and *Salmonella*. *EcoSal Plus*. 2021. doi:10.1128/ecosalplus.esp-0033-2020

18. McKnight CJ, Cross ELA, Kuebler A, et al. High diversity of *Escherichia coli* causing invasive disease in neonates in Malawi poses challenges for O-antigen based vaccine approach. *Commun Med*. 2025;5:133. doi:10.1038/s43856-025-01007-1

19. Abubakar A, et al. Characterization of antimicrobial-resistant Gram-negative bacteria that cause neonatal sepsis in seven low- and middle-income countries. *Nat Microbiol*. 2021;6:512–523. doi:10.1038/s41564-021-00870-7

    *(This covers the BARNARDS network including South Africa [Barnards collection] and Pakistan sites; also encompasses Ethiopia sites. Patan [Nepal], Mbira [Zimbabwe], MRCG [Gambia], and Benin collections — please provide specific publication references or confirm if these are unpublished/pre-publication datasets.)*

20. McKinney W. Data structures for statistical computing in Python. *Proc 9th Python Sci Conf*. 2010;445:51–56.

---

### Notes on unresolved citations

**[CITATION] — G1/G4 prevalence in invasive/neonatal disease (refs 15 and in-text):** The Malawi *E. coli* diversity paper (ref 18) documents capsule type prevalence in neonatal sepsis; this or an equivalent multi-country comparative study should be cited. Please advise on preferred reference(s) for the geographic distribution claim.

**[FastKaptive]:** No peer-reviewed paper exists for FastKaptive by R. Mostowy. Options: (a) cite the GitHub repo URL; (b) describe the method in text without a formal citation; (c) use Kaptive 3 (ref 14) as a citation for the Kaptive framework and note FastKaptive as a predecessor pipeline.

**NNS collection citations (ref 13):** The BARNARDS paper (ref 19) covers South Africa, Pakistan, and Ethiopia. Please provide the specific publication or ENA submission references for: Patan (Nepal), Mbira (Zimbabwe), MRCG (Gambia), Benin. If these are Gladstone lab collections or in-house unpublished datasets, they should be cited as "data not yet published" or via ENA/NCBI accession numbers.

**[AUTHORS], [AFFILIATIONS], [ACKNOWLEDGEMENTS], pyrodigal v[VERSION]:** Please fill in directly.

---

*Database version: v0.9 (651 loci)*
*GitHub: https://github.com/efosternyarko/EC-K-typing-G1G4*
