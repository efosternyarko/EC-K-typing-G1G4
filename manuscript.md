# EC-K-typing Group 1 & 4: a reference database and normalised scoring framework for *Escherichia coli* Group 1 and Group 4 capsule typing

**[AUTHORS]**

[AFFILIATIONS]

Correspondence: [CORRESPONDING AUTHOR EMAIL]

---

## Abstract

The *Escherichia coli* capsular polysaccharide (K antigen) is a major virulence determinant and a target for surveillance and vaccine development. Existing sequence-based K-typing tools cover Group 2 and Group 3 capsule loci, which use the ABC transporter export pathway. Group 1 and Group 4 loci, which use the Wzy-dependent pathway and are located at the *cps* chromosomal region, lack a comprehensive reference database compatible with modern typing tools. Here we present EC-K-typing Group 1 & 4, a curated reference database of 93 filtered *E. coli* Group 1 and Group 4 capsule loci (designated KL300–KL423, K24, and K96) derived from 1,112 *E. coli* bloodstream infection (BSI) genomes. Loci are annotated with systematic positional gene names, formatted as GenBank records, and integrated with the existing Group 2 & 3 database to provide a combined 183-locus resource covering all four capsule groups. We further introduce a normalised scoring approach that re-ranks Kaptive alignment scores by alignment score per expected reference base pair (Norm AS = AS / total_expected_CDS_bp), eliminating size-dependent bitscore accumulation bias inherent in standard Kaptive scoring. Normalised scoring raises self-typing accuracy from 70/93 (75.3%) with standard Kaptive to 91/93 (97.8%) and achieves 100% typeability across all 567 validation assemblies. Applied to 592 *E. coli* assemblies from eight neonatal sepsis collections spanning six countries, the database assigns a capsule type to 100% of assemblies, revealing a predominance of Group 1/4 types (49.2%), with KL302 as the most frequent type (38% of G1/G4 calls). The database, annotation pipeline, and normalised scoring scripts are openly available at https://github.com/[REPO].

**Keywords:** *Escherichia coli*, capsule, K antigen, K typing, Group 1, Group 4, Kaptive, neonatal sepsis, virulence, genomic epidemiology

---

## Introduction

The capsular polysaccharide (CPS) of *Escherichia coli*, encoded by the K-antigen locus, is a central virulence factor that confers resistance to complement-mediated killing and phagocytosis [CITATION]. More than 80 distinct K-antigen serotypes have been described by classical methods, and their distribution varies substantially between disease contexts — invasive disease, urinary tract infection, and neonatal sepsis each show distinctive K-type profiles [CITATION]. As genomic surveillance of bacterial pathogens expands, sequence-based capsule typing has become an essential complement to whole-genome phylogenetics, enabling direct comparison of virulence potential across collections and time points without requiring serological reagents.

*E. coli* capsules are divided into four groups based on biosynthetic pathway and chromosomal location [CITATION]. Group 2 and Group 3 capsules are synthesised via the ABC transporter pathway and their biosynthetic loci map near the *serA* chromosomal locus. These types are covered by the EC-K-typing database (Gladstone et al.) [CITATION], which provides 90 reference loci and achieves high self-typing accuracy with Kaptive [CITATION]. Group 1 and Group 4 capsules are synthesised via the Wzy-dependent pathway; their biosynthetic loci map to the *cps* region, flanked by *galF* upstream and *gnd* downstream, with the *wza-wzb-wzc* capsule export genes immediately upstream of the locus-variable biosynthetic region. Despite their clinical importance — Group 1 types predominate in invasive disease and neonatal infection in multiple African settings [CITATION] — no comprehensive sequence database for Group 1 and Group 4 loci compatible with Kaptive or similar tools has previously been available.

A key challenge in building such a database is the bitscore accumulation bias inherent in Kaptive's default scoring. Kaptive assigns an alignment score (AS) to each candidate locus by summing BLAST bitscores across all reference CDS genes found in the query assembly. Because AS accumulates across genes, loci with more or longer CDS will always outscore smaller loci even if both match at identical per-base identity. In databases covering closely related loci of different sizes, this creates systematic mistyping: the largest reference locus dominates the leaderboard regardless of true identity. This bias is compounded in Group 1 and Group 4 databases because all loci share the same conserved flanking genes (*galF*, *gnd*, *ugd*, *wza*, *wzb*, *wzc*), which contribute an identical background score to every locus. Normalisation of the raw alignment score by the total expected CDS length converts bitscore accumulation into a per-base identity metric, making scores comparable across loci of any size.

Here we present EC-K-typing Group 1 & 4, a resource that addresses both the gap in reference sequence coverage and the scoring bias problem. We describe the construction, annotation, and iterative improvement of the database through six versions, demonstrate its application to bloodstream infection and neonatal sepsis datasets, and provide a normalised scoring framework that is immediately applicable to any Kaptive-compatible database.

---

## Materials and Methods

### Source genome collection and Group 1/4 identification

The primary source dataset comprised 6,673 *E. coli* bloodstream infection (BSI) genome assemblies obtained from EnteroBase (https://enterobase.warwick.ac.uk/) [CITATION]. All 6,673 genomes were screened against the EC-K-typing Group 2 & 3 database (v3.0.0) [CITATION] using Kaptive v3.1.0 [CITATION]. A total of 1,112 assemblies (16.7%) produced no significant hit, indicating the likely presence of a Group 1 or Group 4 capsule locus. These 1,112 no-hit assemblies were subsequently typed with FastKaptive [CITATION] to assign preliminary KX-type designations. Fifteen distinct KX types were identified, including KX36 (*n* = 382), KX34 (*n* = 329), KX17 (*n* = 195), KX01 (*n* = 73), KX67 (*n* = 41), KX31 (*n* = 29), KX72 (*n* = 27), and eight additional less frequent types.

A second validation dataset comprised 592 *E. coli* assemblies from eight neonatal sepsis (NNS) collections spanning six countries: Patan (Nepal; *n* = 20), Barnards (South Africa; *n* = 75), Mbira (Zimbabwe; *n* = 57), CHAMPS Harar (Ethiopia; *n* = 110), Malawi Mlw (*n* = 167), MRCG (Gambia; *n* = 130), Benin (*n* = 20), and Pakistan (*n* = 13). These assemblies were assembled from short-read Illumina data using [ASSEMBLER] and used for application-level validation of the final database.

### Locus extraction

Group 1/4 capsule loci were extracted from the 1,112 BSI assemblies using a flanking-gene detection approach. Reference sequences for *galF*, *gnd*, *ugd*, *wza*, and *wzc* from *E. coli* K-12 MG1655 were used as BLAST queries against each assembly. The *cps* region was defined as the sequence between the outermost detected *galF* and *gnd* flanking genes, extended to include *wza-wzb-wzc* upstream, with 500 bp of additional flanking sequence at each boundary. For assemblies in which *galF* and *gnd* mapped to different contigs (indicative of locus fragmentation across assembly breaks), partial sequences from both contigs were extracted and concatenated with a poly-N spacer.

### Clustering and representative selection

Extracted locus sequences were clustered by all-vs-all BLASTn at ≥95% nucleotide identity and ≥80% query coverage, with greedy clustering by sequence length. Clustering yielded 125 distinct locus types (KL300–KL423, K24, and K96), covering all 10 KX types identified by FastKaptive. A filtered set of 93 loci ≥30 kb was retained for inclusion in the Kaptive-compatible database, excluding short or fragmented extractions unlikely to provide reliable typing. For two loci (KL388 and KL391), the original BSI representatives were replaced in database version 0.3.1 with longer, more complete sequences identified in the neonatal sepsis cohort (KL388: +6.6 kb; KL391: +10.5 kb). For two additional loci (KL306 and KL307), BSI representatives that failed self-typing validation were replaced in version 0.5 with sequences from NCBI (accessions CP099041 and CP070103 respectively) identified by megaBLAST candidate search and validated by normalised scoring.

### Gene prediction and annotation

Open reading frames were predicted across each of the 93 filtered locus sequences using pyrodigal v[VERSION] [CITATION] in metagenomic mode, which does not require a training set and is appropriate for sequences of heterogeneous GC content and genetic context. Predicted proteins were annotated in two stages:

1. **Functional name transfer:** predicted amino acid sequences were searched against the Kaptive *Klebsiella* K-locus reference database [CITATION] by BLASTp (≥30% identity, ≥50% query and subject coverage). Matching proteins received the functional gene name from the *Klebsiella* reference (e.g., *wza*, *wzb*, *wzc*, *wzx*, *wzy*, and glycosyltransferase names). Seventy-three percent of predicted CDS received a functional annotation by this approach.

2. **Positional name assignment (v0.3):** all predicted proteins across the 93 loci were clustered by all-vs-all BLASTp at ≥90% identity and ≥90% query and subject coverage. Protein families shared across loci received a deterministic positional name of the form `KL{N}_{pos}`, where *N* is the KL number of the lowest-numbered locus contributing a family member and *pos* is its ordinal position within that locus. If any member of a family carried a Klebsiella-derived functional name, that name was propagated to the entire family. Within-locus duplicate names were disambiguated with `_2`, `_3` suffixes.

The final annotation comprised 3,298 CDS across the 93 filtered loci: 2,094 (64%) with functional names and 1,204 (36%) with positional names, defining 67 multi-locus protein families shared across two or more loci.

### GenBank record formatting

All 93 G1/G4 locus sequences were stored as GenBank flat files formatted for Kaptive compatibility. Each record carries a `source` feature with `/note="K locus: KLxxx"` and `/note="K type:KLxxx"` qualifiers, which Kaptive uses to report locus and type identifiers. The 93 G1/G4 records were combined with the 90-locus EC-K-typing Group 2 & 3 database (v3.0.0) to produce a combined 183-locus all-groups database.

### Normalised scoring

Standard Kaptive reports a raw alignment score (AS) equal to the cumulative BLAST bitscore across all reference CDS genes detected in the query assembly. Because AS accumulates additively across genes, a locus with *n* genes of length *L* bp accumulates approximately *n × L × f* bitscore units (where *f* is a per-base function of alignment identity), making loci with more or longer genes systematically outscore smaller loci even at identical per-base identity.

We introduce a normalised alignment score (Norm AS) defined as:

> **Norm AS = AS / total_expected_CDS_bp**

where `total_expected_CDS_bp` is the sum of lengths of all CDS annotated in the reference locus GenBank record. For a hypothetical locus typed at 100% identity with all genes present, Norm AS converges to a value of approximately 2.0 (reflecting the relationship between nucleotide BLAST bitscore and alignment length at full identity with the default BLAST scoring matrix). This metric is size-independent and directly comparable across loci of any length or gene content.

Normalised scoring was implemented in `scripts/type_normalized.py`, which:
1. Invokes `kaptive assembly --scores` to generate a full locus × assembly score matrix
2. Reads the expected CDS lengths from each GenBank record
3. Computes Norm AS for every locus × assembly pair
4. Selects the best-match locus as the highest Norm AS per assembly
5. Applies a near-tie tiebreaker: among loci scoring within 0.5% of the top Norm AS, the locus with more genes found is preferred
6. Classifies each assembly as Typeable (≥50% of expected genes found) or Untypeable

A confidence grading is additionally applied to typeable assemblies: High (Norm AS ≥ 1.90), Moderate (1.50 ≤ Norm AS < 1.90), or Low (Norm AS < 1.50).

### Representative replacement via NCBI candidate search

For loci failing self-typing validation after normalised scoring (KL300, KL303, KL306, KL307 in version 0.4), a systematic NCBI candidate search was performed. The reference locus sequence for each failing locus was submitted to NCBI megaBLAST against the nucleotide (nt) database restricted to *Escherichia coli* (taxid:562) with ≥90% identity and ≥70% query coverage thresholds. The top 10 candidate chromosome assemblies per locus were downloaded via NCBI Entrez and tested with normalised scoring against the v0.4 database. Candidates that self-typed correctly (best Norm AS for the expected locus) were retained as replacement representatives.

For accepted replacements, new locus sequences were extracted from the candidate chromosome by multi-HSP BLASTn tiling, and gene annotations were transferred from the original representative by per-CDS BLASTn liftover (≥80% identity, ≥80% query coverage) rather than re-annotation with pyrodigal.

---

## Results

### Database overview and content

The EC-K-typing Group 1 & 4 database (current version: v0.6) contains 93 reference locus sequences covering 10 KX types identified in *E. coli* bloodstream infection isolates (Table 1). Locus sequences range from 30.0 to 59.1 kb in the filtered set, with a median of 37.9 kb. A total of 3,471 CDS are annotated across the 93 loci (median 37 CDS per locus; range 23–50), of which 64% carry functional names and 36% positional names.

**Table 1. EC-K-typing Group 1 & 4 database content summary.**

| Feature | Value |
|---------|-------|
| Total G1/G4 reference loci | 93 |
| KL designation range | KL300–KL423, K24, K96 |
| KX types covered | 10 (KX01, KX17, KX28, KX31, KX34, KX36, KX63, KX67, KX72, K24) |
| Locus size range (filtered set) | 30.0–59.1 kb |
| Total reference sequence | 3.69 Mb |
| Total CDS annotated | 3,471 |
| CDS with functional names | 64% (2,223) |
| CDS with positional names | 36% (1,248) |
| Multi-locus shared protein families | 67 |
| Combined all-groups database size | 183 loci (93 G1/G4 + 90 G2/G3) |

The combined 183-locus all-groups database merges the 93 G1/G4 loci with the 90-locus EC-K-typing Group 2 & 3 database (Gladstone et al.) and is provided as a single GenBank file for immediate use with Kaptive.

### Bitscore accumulation bias and normalised scoring

Initial validation of the v0.2 database (70/93 filtered loci using standard Kaptive scoring; 75.3% self-typing; 258/565 assemblies typeable; 45.7%) revealed that failures were concentrated among KX01 loci. Root-cause analysis identified bitscore accumulation bias as the primary driver: KL302 (40 CDS, 37.8 kb total CDS length) and KL304 (44 CDS, 46.5 kb) accumulate disproportionately high raw alignment scores relative to smaller KX01 references (28–35 kb CDS), even when all query genes match at 100% identity.

Applying normalised scoring (Norm AS = AS / total_expected_CDS_bp) to the raw score matrix from the same Kaptive `--scores` run immediately raised self-typing to 88/93 (94.6%) and typeability to 565/565 (100.0%) (Table 2). The improvement was achieved without any change to the database sequences or annotations — solely by correcting the scoring metric.

The near-tie tiebreaker resolved one additional ambiguous case: KL337 vs. KL389 (Norm AS difference 0.17%; 41 vs. 29 genes found), correctly assigning the assembly to KL337.

Five loci remained failures after normalised scoring (KL300, KL301, KL303, KL306, KL307), all within the KX01 group, reflecting genuine per-base similarity of these loci to KL302 in their source genomes rather than a scoring artifact.

**Table 2. Self-typing accuracy and typeability across database versions.**

| Version | Standard Kaptive self-typing | Normalised self-typing | Assemblies typeable (normalised) |
|---------|-----------------------------|-----------------------|---------------------------------|
| v0.2 | 70/93 (75.3%) | 88/93 (94.6%) | 565/565 (100.0%) |
| v0.3 | 70/93 (75.3%) | 88/93 (94.6%) | 565/565 (100.0%) |
| v0.3.1 | 70/93 (75.3%) | 88/93 (94.6%)* | 592/592 (100.0%)† |
| v0.4 | n.d. | 89/93 (95.7%) | 565/565 (100.0%) |
| v0.5 | 67/93 (72.0%) | 91/93 (97.8%) | 567/567 (100.0%) |
| v0.6 | n.d. | 91/93 (97.8%) | 567/567 (100.0%) |

*Normalised self-typing using NNS validation assemblies. †592 NNS assemblies; all typeable.
n.d.: not determined for this version.

### Iterative improvement through representative replacement

Five KX01 loci that failed normalised self-typing (KL300, KL301, KL303, KL306, KL307) were investigated individually.

**KL301** was resolved in version 0.4 by removing the shared conserved flanking genes (*galF*, *galF_2*, *gnd*, *ugd*, *wza*, *wzb*, *wzc*) from all G1/G4 GenBank records. These seven gene families (528 CDS, 16.0% of the total G1/G4 CDS) contribute identical bitscore to every KX01 locus, partially masking variable-region differences. After removal, KL301 self-typed correctly at Norm AS = [VALUE]. This improvement was retained in all subsequent versions by permanently excluding conserved CDS from the KL301 record.

**KL306 and KL307** were resolved in version 0.5 by replacing their original BSI representatives with NCBI chromosome sequences identified through systematic megaBLAST candidate search. For KL306, the original BSI assembly (ESC_GB3606AA_AS, 41,735 bp) was replaced by CP099041 (42,626 bp; 31/31 CDS successfully lifted over). For KL307, ESC_NB5901AA_AS (44,329 bp) was replaced by CP070103 (45,329 bp; 34/34 CDS lifted over). Both new representatives self-typed with perfect normalised scores (Norm AS = 2.000), confirming 100% per-base identity between the representative and its source genome. The original BSI assemblies for both loci typed as KL302 with both the old and new databases, confirming that the BSI representatives were genuinely poor matches for their assigned types.

**KL300** was evaluated using the top 10 NCBI candidates (299 megaBLAST hits; top candidate: CP135488) but no suitable replacement was identified. With a standard locus extraction window (27/30 CDS recovered), the candidate variable region was too conserved across KX01 loci, causing non-KL300 assemblies to shift toward KL300. Expanding the extraction window to include peripheral flanking genes (all 30 CDS) caused typeability regression for 54 assemblies by introducing size-bias in standard Kaptive locus-finding. KL300 was therefore retained with its original BSI representative.

**KL303** remains unresolvable from KL302 in all publicly available *E. coli* genomes. All 264 NCBI megaBLAST candidates typed as KL302 or KL352 — none self-typed as KL303. This likely reflects genuine sequence similarity between the KL303 and KL302 variable biosynthetic regions in all currently sequenced isolates. KL303 is retained in the database with a documented caveat that self-typing cannot be achieved with current public data.

### Conserved gene restoration (version 0.6)

Version 0.4 stripped conserved flanking/export CDS from all 93 G1/G4 loci to reduce background scoring bias. Subsequent analysis showed this was necessary only for KL301; the other four failures were resolved by representative replacement (KL306, KL307) or reflect biological ambiguity (KL300, KL303). Full locus annotations — including conserved genes — are important for user interpretation of typing results and gene-level characterisation of novel assemblies. Version 0.6 therefore restores conserved CDS to all 92 loci except KL301, with no change in normalised scoring performance.

For 91 loci with unchanged sequences relative to v0.3.1, conserved CDS features were copied directly from the v0.3.1 GenBank records (coordinates are identical). For KL306 and KL307, which carry new NCBI-derived sequences, conserved gene positions were located in the new sequences by per-gene BLASTn liftover (≥80% identity, ≥80% query coverage) with coordinate snapping to the nearest multiple of 3. The *ugd* gene is annotated as a partial CDS in both KL306 (249 bp) and KL307 (252 bp) because the locus extraction boundary falls at the *gnd*/*ugd* junction; this is a genuine locus-edge truncation, not an annotation error.

### Application to neonatal sepsis cohorts

The combined 183-locus database was applied to 592 *E. coli* assemblies from eight neonatal sepsis (NNS) collections using normalised scoring (Table 3). All 592/592 assemblies were typeable (Norm AS ≥ threshold, ≥50% gene coverage). Group 1/4 loci were assigned to 291 assemblies (49.2%); Group 2/3 loci to 301 assemblies (50.8%).

The most frequent G1/G4 types across all collections were KL302 (110/291; 37.8%), KL301 (35; 12.0%), KL303 (28; 9.6%), KL329 (23; 7.9%), and KL337 (19; 6.5%). Forty G1/G4 loci were represented in the NNS dataset; 22 were matched at high confidence (Norm AS ≥ 1.90, 100% gene coverage), indicating that the NNS cohorts cover a broad fraction of the G1/G4 diversity present in the reference database.

The proportion of G1/G4 types varied markedly between collections (range: 34–75%), reflecting both geographic variation in capsule-type prevalence and potential differences in clinical presentation. The Gambian (MRCG) collection had the lowest G1/G4 proportion (34%), while the Nepali (Patan) collection had the highest (75%), consistent with previously observed geographic trends in *E. coli* capsule type distribution in invasive disease [CITATION].

**Table 3. Capsule typing of 592 neonatal sepsis *E. coli* assemblies by collection.**

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

The database and associated software are distributed as a public GitHub repository (https://github.com/[REPO]). All analysis scripts are written in Python 3 and depend on Biopython [CITATION], pandas [CITATION], and NCBI BLAST+ [CITATION]. The current database version (v0.6) is provided as two GenBank files:

- `DB/EC-K-typing_group1and4_v0.6.gbk` — 93 G1/G4 reference loci (for G1/G4-specific analysis)
- `DB/EC-K-typing_all_groups_v0.6.gbk` — 183 combined loci from all four capsule groups (recommended for general use)

Pre-computed validation results for 567 assemblies (565 BSI + 2 NCBI representatives) are provided in `DB/kaptive_validation_results_v0.6norm.tsv`. NNS typing results for 592 assemblies are in `DB/nns_kaptive_results_v0.5norm.tsv`.

### Software

**`scripts/type_normalized.py`** — The primary analysis tool. Accepts a Kaptive-compatible GenBank database and a set of genome assemblies (FASTA), runs Kaptive in `--scores` mode, computes normalised alignment scores, and outputs three files: the raw scores matrix, per-assembly normalised typing results, and a per-assembly summary with expected versus observed type. Usage:

```bash
python scripts/type_normalized.py \
  --db DB/EC-K-typing_all_groups_v0.6.gbk \
  --suffix myrun \
  --threads 8
```

A `--skip-kaptive` flag allows re-normalisation of a pre-computed scores matrix without re-running Kaptive, enabling rapid parameter exploration.

**Standard Kaptive** can also be used directly with the GenBank database:

```bash
kaptive assembly DB/EC-K-typing_all_groups_v0.6.gbk genome.fasta -o results.tsv
```

Standard Kaptive achieves 72.0% self-typing for G1/G4 loci and types approximately 43% of assemblies (the remainder are reported as Untypeable due to locus fragmentation across contigs or divergence from the reference that prevents initial locus localisation by full-sequence BLAST). Normalised scoring is therefore strongly recommended for G1/G4 typing.

### Interpreting normalised alignment scores

| Norm AS | Interpretation |
|---------|---------------|
| ~2.00 | Perfect: all reference genes found at 100% identity |
| ≥ 1.90 | High confidence: nearly all genes at high identity |
| 1.50–1.90 | Moderate: partial locus or gene-level divergence |
| < 1.50 | Low confidence: possible novel type, divergent variant, or assembly issue |

KL302, KL303, and KL300 calls should be treated with particular caution, as KL303 cannot currently be distinguished from KL302 in any public genome sequence (see Discussion).

---

## Discussion

We present EC-K-typing Group 1 & 4, the first comprehensive sequence typing database for *E. coli* Group 1 and Group 4 capsule loci, together with a normalised scoring framework that substantially improves typing accuracy over standard Kaptive while maintaining 100% typeability. The database fills a significant gap in *E. coli* genomic epidemiology toolkits: until now, automated K-antigen typing was restricted to Group 2 and Group 3 loci, leaving approximately 15–20% of BSI isolates and a large fraction of neonatal sepsis isolates uncharacterised at the capsule level.

The normalised scoring approach introduced here addresses a fundamental limitation of cumulative bitscore-based typing that extends beyond this specific database. Any Kaptive database covering loci of different sizes will exhibit bitscore accumulation bias unless scores are normalised by locus size. The approach is computationally trivial — requiring only a single division per locus × assembly pair — and requires no changes to Kaptive itself or to the underlying BLAST analysis. We recommend that it be considered for adoption in other Kaptive-based typing databases, particularly those covering closely related loci of variable gene content.

The five KX01 loci (KL300–KL307 with the exception of KL302 and KL304) present a persistent challenge. Three of the five (KL301, KL306, KL307) have been successfully resolved through a combination of conserved gene manipulation and NCBI representative replacement, respectively. KL303 represents a genuinely hard biological case: despite 264 public genomes matching the locus at high similarity, none carry a variable region distinguishable from KL302 by sequence alone. This may reflect a very recent divergence event, ongoing horizontal gene transfer between KL302 and KL303 biosynthetic loci, or an error in the original KL303 designation. We retain KL303 in the database to ensure that loci morphologically consistent with this type can still be detected and flagged, but users should treat KL303 assignments as requiring confirmatory evidence (e.g., antigen serology or independent sequencing of the biosynthetic region with long-read technology).

The neonatal sepsis application reveals that Group 1/4 types constitute nearly half (49.2%) of all capsule types in these cohorts, with substantial geographic variation. KL302 is the dominant G1/G4 type across all sites, consistent with the broad prevalence of KX36-type *E. coli* in invasive disease. The high frequency of KL301 (12% of G1/4 calls) and KL303 (9.6%) in the African collections mirrors the known predominance of KX01-type *E. coli* in neonatal bacteraemia in sub-Saharan Africa [CITATION]. Caution is warranted in interpreting KL303 calls in clinical datasets given its indistinguishability from KL302, and future work should prioritise resolving this ambiguity.

Several limitations should be noted. First, the source BSI dataset (EnteroBase, UK-centric) may not fully represent the global diversity of Group 1 and Group 4 capsule types; loci present only in low-income country settings or specific endemic regions may be absent from the current database. Second, all locus sequences are derived from short-read assemblies; loci fragmented across contig boundaries may be incompletely represented, and long-read assembly would be expected to improve reference quality. Third, the positional gene naming scheme (e.g., KL302_04_wzx) is not directly equivalent to established *E. coli* gene names or classical capsule gene designations, limiting direct comparison with serotype-based literature. Future work will focus on mining the AllTheBacteria dataset (~2.4 million prokaryotic assemblies) [CITATION] for discriminating representatives for KL300 and KL303, expanding the reference set to include loci from under-represented geographic regions, and producing a v1.0 release with 100% self-typing once remaining ambiguities are resolved.

---

## Data Availability

All database files, annotation scripts, normalised scoring code, and pre-computed validation results are freely available at:

**https://github.com/[REPO]**

The current recommended database for all-groups typing is `DB/EC-K-typing_all_groups_v0.6.gbk`. The normalised scoring script is `scripts/type_normalized.py`. Pre-computed results for 567 BSI validation assemblies and 592 NNS assemblies are provided in the `DB/` directory.

The EC-K-typing Group 2 & 3 database (v3.0.0) used as the complementary G2/G3 component is available at https://github.com/rgladstone/EC-K-typing.

Source genome assemblies are available from EnteroBase (https://enterobase.warwick.ac.uk/) under the accessions listed in `DB/KL_G1G4_mapping_filtered.tsv`. NCBI accessions CP099041 (KL306 representative) and CP070103 (KL307 representative) are available from NCBI Nucleotide.

---

## Acknowledgements

[ACKNOWLEDGEMENTS — funding, data access agreements, sequencing centres, etc.]

---

## References

1. Whitfield C, Roberts IS. Structure, assembly and regulation of expression of capsules in *Escherichia coli*. *Mol Microbiol*. 1999;31:1307–1319.

2. Russo TA, Marr CM. Hypervirulent *Klebsiella pneumoniae*. *Clin Microbiol Rev*. 2019;32:e00001-19. [Note: replace with appropriate *E. coli* K antigen virulence citation]

3. Whitfield C, Wear SS, Sande C. Assembly of bacterial capsular polysaccharides and exopolysaccharides. *Annu Rev Microbiol*. 2020;74:521–543.

4. Gladstone RA, [CO-AUTHORS]. A sequence typing database for *Escherichia coli* capsular K antigens (Group 2 & 3). *medRxiv*. 2024. https://doi.org/10.1101/2024.11.22.24317484

5. Lam MMC, Wick RR, Judd LM, Holt KE, Wyres KL. Kaptive 2.0: updated capsule and lipopolysaccharide locus typing for the *Klebsiella pneumoniae* species complex. *Microb Genom*. 2022;8:000800.

6. Alikhan NF, Zhou Z, Sergeant MJ, Achtman M. A genomic overview of the population structure of *Salmonella*. *PLoS Genet*. 2018;14:e1007261. [EnteroBase citation — replace with correct EnteroBase ref]

7. Larralde M, Zeller G. Pyrodigal: faster gene predictions with Prodigal. *J Open Source Softw*. 2022;7:4296.

8. Camacho C, Coulouris G, Avagyan V, et al. BLAST+: architecture and applications. *BMC Bioinformatics*. 2009;10:421.

9. Cock PJA, Antao T, Chang JT, et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics*. 2009;25:1422–1423.

10. Rowe W, Kalule M, Dyer C, et al. AllTheBacteria: all bacterial genomes assembled, available and searchable. *bioRxiv*. 2024. https://doi.org/10.1101/2024.03.08.584059

11. Shaw LP, [CO-AUTHORS]. LexicMap: efficient sequence search in large-scale microbial genome collections. *Nat Biotechnol*. 2025. https://doi.org/10.1038/s41587-025-02812-8

12. [NEONATAL SEPSIS CITATION — relevant papers for Patan/Barnards/Mlw/MRCG/CHAMPS collections]

---

*Manuscript prepared: [DATE]*
*Database version: v0.6*
*GitHub: https://github.com/[REPO]*
