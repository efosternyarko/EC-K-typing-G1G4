# EC K-typing G1/G4 Database Build Protocol

**Repository:** https://github.com/efosternyarko/EC-K-typing-G1G4
**Author:** Dr Ebenezer Foster-Nyarko
**Current version:** v0.9 (651 loci; 651/651 self-typing 100%)
**Last updated:** March 2026

---

## Overview

This protocol describes the construction and validation of a Kaptive-compatible
reference database for typing *Escherichia coli* Group 1 and Group 4 (G1/G4)
capsular polysaccharide (K) loci. G1/G4 loci use the Wzy-dependent pathway
and are flanked by the conserved chromosomal genes *galF* (upstream) and *gnd*
(downstream).

The database is intended for use **after** screening with the G2/G3 database
(sequential typing): run G2/G3 first; apply this G1/G4 database only to
untypeable assemblies.

> **Why sequential, not combined?**
> A combined all-groups database causes wzy-interference: G1/G4 (wzy-based)
> loci consistently outscore G2/G3 (kps-based) loci in direct competition,
> rendering genuine G2/G3 loci untypeable. Validation against the Malawi
> dataset confirmed this — 82/100 correctly G2/G3-typed assemblies became
> untypeable when the combined database was used.

---

## Prerequisites

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| Kaptive | ≥ 3.2.1 | K-locus typing |
| minimap2 | ≥ 2.28 | Sequence alignment |
| BLAST+ | ≥ 2.15.0 | Flanking gene screening, BLASTp annotation |
| MMseqs2 | 17-b804f | Clustering novel loci |
| pyrodigal | ≥ 3.0 | Gene prediction |
| Biopython | ≥ 1.81 | GenBank parsing/writing |
| LexicMap | 0.8.1 | Targeted genome search (AllTheBacteria) |

Install Python dependencies:
```bash
pip install kaptive pyrodigal biopython pandas
```

### Data resources

- **AllTheBacteria (ATB) v202408** — 888 batch archives at
  `/scratch2/js66/atb_archives/atb.assembly.incr_release.202408.batch.NNN.tar.xz`
  on M3 HPC (Monash). Each archive contains FASTA assemblies named
  `<BioSample_accession>.fa`.
- **ATB LexicMap index** — `s3://allthebacteria-lexicmap/202408/` (public,
  AWS eu-west-2). Searched via LexicMap on an EC2 instance.
- **Klebsiella K-locus primary reference** — bundled with Kaptive at
  `kaptive/data/Klebsiella_k_locus_primary_reference.gbk`. Used for BLASTp
  gene name transfer.

### Cluster / cloud access

- **M3 HPC (Monash):** `ebenezef@m3.massive.org.au`
  Scratch space: `/home/ebenezef/js66_scratch/ebenn/`
- **AWS:** eu-west-2 region, credentials configured via `aws configure`

### Repository layout

```
EC-K-typing-G1G4/
├── DB/                          # Reference databases and validation outputs
│   ├── EC-K-typing_group1and4_v0.8.gbk
│   ├── reference_loci_v0.6.fasta
│   ├── kaptive_validation_summary_v0.8norm.tsv
│   └── atb_lexicmap_results/
├── flanking_genes/
│   └── flanking_genes.fasta     # galF + gnd reference sequences
└── scripts/                     # All pipeline scripts (described below)
```

---

## Phase 1 — ATB Screening for G1/G4 Candidates

**Objective:** Identify ATB assemblies containing G1/G4 *cps* loci by
detecting the conserved flanking genes *galF* and *gnd*.

**HPC location:**
```
/home/ebenezef/js66_scratch/ebenn/atb_screen/
├── results/          # Per-batch BLASTn hit files (batch_N_hits.tsv)
├── novel_loci/       # Extracted novel locus FASTAs
├── logs/
├── atb_g1g4_candidates.tsv
└── batches_with_candidates.txt
```

### Step 1.1 — Screen all 888 ATB batches (SLURM array)

**Script:** `scripts/submit_atb_screen.sh` → calls `scripts/atb_screen_batch.py`

```bash
# On M3:
cd /home/ebenezef/js66_scratch/ebenn/atb_screen
sbatch /home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/scripts/submit_atb_screen.sh
```

Each array task (`1–888`, max 50 concurrent) processes one ATB batch:
1. Extracts all assemblies from the `.tar.xz` archive into a temp directory
2. Builds a combined FASTA
3. BLASTn against `flanking_genes/flanking_genes.fasta` (*galF* + *gnd*)
4. Writes `results/batch_N_hits.tsv`

**Key thresholds:** ≥80% nucleotide identity, ≥70% query coverage, e-value ≤1e-10.

**Monitor:**
```bash
squeue -u ebenezef
# Count completed:
ls results/batch_*_hits.tsv | wc -l   # expect 665 of 888 with candidates
```

### Step 1.2 — Combine screening results

**Script:** `scripts/combine_atb_hits.py`

```bash
# On M3:
python3 /home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/scripts/combine_atb_hits.py
```

Filters for genomes with **both** *galF* AND *gnd* hits on the **same contig**
(guarantees the full locus is intact and extractable). Writes:
- `atb_g1g4_candidates.tsv` — genomes to extract loci from
- `batches_with_candidates.txt` — list of batch numbers containing ≥1 candidate
- `atb_screen_summary.txt`

---

## Phase 2 — Novel Locus Extraction

**Objective:** For each G1/G4 candidate genome, extract the locus sequence
between *galF* and *gnd*, discarding sequences that match the existing
reference database (known types).

### Step 2.1 — Build novelty-screening BLAST database

**Script:** `scripts/setup_novelty_db.sh`

```bash
# On M3, run once:
bash /home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/scripts/setup_novelty_db.sh
```

Builds a BLAST nucleotide database from `DB/reference_loci_v0.6.fasta`
(the 93 existing reference loci) at:
`/home/ebenezef/js66_scratch/ebenn/atb_screen/novelty_db/g1g4_refs`

### Step 2.2 — Extract novel loci (SLURM array)

**Script:** `scripts/submit_extract_novel.sh` → calls `scripts/extract_novel_loci.py`

```bash
sbatch /home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/scripts/submit_extract_novel.sh
```

Each array task processes one candidate batch (`1–268`, max 50 concurrent):
1. Loads same-contig candidates for the batch
2. BLASTn each candidate contig against the novelty DB:
   - If ≥95% identity AND ≥80% coverage → **known type, skip**
   - Otherwise → **novel**
3. BLASTn *galF* and *gnd* against the contig to get precise coordinates
4. Extracts the region between *galF* end and *gnd* start, ±500 bp flanking
5. Discards extractions < 15,000 bp (likely partial)
6. Writes `novel_loci/batch_N_novel_loci.fasta`

**Result:** 347,524 novel locus sequences across 268 batches.

---

## Phase 3 — Clustering Novel Loci

**Objective:** Collapse 347,524 raw novel locus sequences into unique KL types
by clustering at 95% identity / 80% bidirectional coverage.

### Step 3.1 — Run MMseqs2 clustering

**Script:** `scripts/cluster_novel_loci.sh`

```bash
sbatch /home/ebenezef/js66_scratch/ebenn/EC-K-typing-G1G4/scripts/cluster_novel_loci.sh
```

MMseqs2 parameters:
```
--min-seq-id 0.95    # 95% nucleotide identity
-c 0.80              # 80% coverage
--cov-mode 0         # bidirectional coverage
--cluster-mode 0     # greedy set cover (default)
--threads 16
```

### Step 3.2 — Assign KL type numbers

The clustering script embeds a Python block that assigns KL numbers to
cluster representatives, sorted by descending cluster size (most common
types get the lowest numbers).

> **Important:** The script header comments mention `KL344` as the start,
> but the existing database (v0.6) already contains loci up to **KL391**.
> The correct start number is **KL392**. This was corrected by running an
> updated `assign_novel_kl_types.py` script on M3:

```python
# assign_novel_kl_types.py (run separately on M3 after clustering)
KL_START = 392   # correct start — v0.6 highest is KL391
```

```bash
# On M3 (after pip install biopython --user):
python3 ~/assign_novel_kl_types.py
```

**Result:** 577 novel KL types assigned (KL392–KL968).
Output files:
- `atb_clustering/novel_kl_types.fasta` — representative sequences
- `atb_clustering/novel_kl_summary.tsv` — KL type, accession, size, cluster count

---

## Phase 4 — Annotation and Database Build

**Objective:** Predict and annotate CDS in novel loci, then merge with the
existing v0.6 database to produce v0.8.

### Step 4.1 — Run annotation and merge

**Script:** `scripts/annotate_and_merge_atb_loci.py`

```bash
# Local (macOS), requires pyrodigal, biopython, blastp:
python3 scripts/annotate_and_merge_atb_loci.py 2>&1 | tee /tmp/annotate_v08.log
```

Pipeline steps within the script:
1. **Partial-capture screening:** BLASTn each novel locus against the existing
   v0.6 reference. Reject if ≥80% of the novel sequence aligns at ≥95%
   identity to any existing locus (catches ATB assemblies that are partial
   captures of existing reference sequences).
2. **CDS prediction:** pyrodigal in metagenomic mode (`meta=True`)
3. **BLASTp annotation:** predicted proteins vs Klebsiella K-locus primary
   reference; transfer gene names at ≥30% identity, ≥50% query coverage
4. **GenBank record construction:** Kaptive-compatible format with `source`
   feature containing `note: "K locus: KLxxx"`
5. **Merge:** appends annotated novel records to v0.6 GenBank → v0.8

**Key paths:**
```python
NOVEL_FASTA  = "atb_clustering/novel_kl_types.fasta"
EXISTING_GBK = "DB/EC-K-typing_group1and4_v0.6.gbk"
OUTPUT_GBK   = "DB/EC-K-typing_group1and4_v0.8.gbk"
BLASTP       = "/Users/lshef4/lshef4_to_copy/anaconda3/envs/my_python.env/bin/blastp"
KLEB_REF_GBK = "/Users/lshef4/Library/Python/3.9/lib/python/site-packages/kaptive/data/Klebsiella_k_locus_primary_reference.gbk"
```

**Result:** `DB/EC-K-typing_group1and4_v0.8.gbk` — 670 records parsed
(93 existing + 577 novel), 57.2% CDS annotated, all passed validation.

> **Note:** Three loci (KL486, KL562, KL812) were subsequently removed
> after self-typing validation revealed they were partial captures of
> KL312, KL391, and KL388 respectively (see Phase 5). Final v0.8: 667 loci.

---

## Phase 5 — Self-typing Validation (Existing 93 Loci)

**Objective:** Confirm that the representative assembly for each of the
original 93 BSI loci correctly types to its expected KL when screened
against the complete v0.8 database.

### Step 5.1 — Run normalised scoring validation

**Script:** `scripts/type_normalized.py`

```bash
python3 scripts/type_normalized.py 2>&1 | tee /tmp/validate_v08c.log
```

This script:
1. Runs `kaptive assembly` with `--scores` against v0.8 on all 93 representative
   FASTAs
2. Calls `normalise_novel_scores.py` to re-rank by `AS / total_expected_gene_bp`
   (normalised score), correcting for locus size bias
3. Outputs a per-assembly validation summary

**Normalised scoring rationale:** Raw Kaptive alignment scores (AS) are
proportional to locus size. A large locus will always outscore a small one
even at lower fractional coverage. Dividing AS by the total expected CDS
length for each locus normalises for size and gives a fair comparison.

**Initial result (first run):** 88/93 correctly typed — failures:
KL300, KL303, KL312, KL388, KL391

### Step 5.2 — Diagnose and fix regressions (KL312, KL388, KL391)

BLASTn analysis revealed:
| Novel locus | Identity to existing | Match |
|-------------|----------------------|-------|
| KL486 | 100% identity, 97% coverage | of KL312 |
| KL562 | 99.8% identity, 100% coverage | of KL391 |
| KL812 | 99.2% identity, 100% coverage | of KL388 |

These three novel loci were partial ATB captures that slipped through the
bidirectional clustering filter (they matched the existing loci but were
shorter, so the coverage check was not symmetric). They were removed from
the v0.8 GenBank:

```bash
# Remove KL486, KL562, KL812 records from v0.8 GenBank, re-run annotation script
# (partial-capture screening step in annotate_and_merge_atb_loci.py now catches these)
```

After removal and re-validation: **91/93 correctly typed, 100% typeability**.
Remaining failures: KL300 and KL303 (see Phase 7).

> **Note on appending:** Kaptive `--scores` appends to an existing TSV by
> default. Always delete the scores file before re-running:
> ```bash
> rm -f DB/kaptive_scores_v0.8norm.tsv
> ```

---

## Phase 6 — Novel Loci Self-typing Validation (M3)

**Objective:** Validate that each of the 574 novel loci (KL392–KL968) can
be correctly typed from its representative ATB assembly using the complete
v0.8 database.

### Step 6.1 — Upload required files to M3

```bash
scp DB/EC-K-typing_group1and4_v0.8.gbk \
    ebenezef@m3.massive.org.au:~/EC-K-typing-G1G4/DB/

scp atb_clustering/novel_kl_summary.tsv \
    ebenezef@m3.massive.org.au:~/EC-K-typing-G1G4/DB/

scp scripts/validate_novel_loci_m3.sh \
    scripts/normalise_novel_scores.py \
    ebenezef@m3.massive.org.au:~/EC-K-typing-G1G4/scripts/
```

### Step 6.2 — Submit validation job

**Script:** `scripts/validate_novel_loci_m3.sh`

```bash
# On M3:
sbatch ~/EC-K-typing-G1G4/scripts/validate_novel_loci_m3.sh
```

The script:
1. Installs Kaptive v3 (`pip install --user "kaptive>=3.0"`)
2. Loads `module load minimap2/2.28`
3. Scans all 268 candidate batch archives (`tar tf`) to locate the 574
   representative accessions. This scanning approach is required because
   cluster representatives are not guaranteed to appear in the
   `batch_*_hits.tsv` files (they may have been identified as novel in
   a different batch's cross-screening step).
4. Extracts FASTAs grouped by archive for efficiency
5. Runs `kaptive assembly --scores` against v0.8
6. Runs `normalise_novel_scores.py` to check self-typing correctness

**Monitor:**
```bash
squeue -u ebenezef
tail -f ~/EC-K-typing-G1G4/logs/val_novel_JOBID.log
```

**Expected result:** ~574/574 correctly self-typed (pending completion).

---

## Phase 7 — Resolving KL300 and KL303 via LexicMap

**Background:** KL300 and KL303 fail self-typing in v0.8. Both loci share the
KX01 polysaccharide gene cluster with KL302, and the original BSI representative
assemblies type to KL302 under normalised scoring. A combined all-groups DB
search cannot resolve these because they type to KL302 in the existing dataset —
ATB genomes genuinely carrying KL300/KL303 would have been excluded from the
novel loci set (they matched KL302 at >95%/80%). A targeted sequence search
across all of ATB is required.

### Step 7.1 — Prepare query FASTA

`DB/kl_failing_query.fasta` contains the existing KL300 (45,380 bp) and
KL303 (40,992 bp) locus sequences extracted from v0.8.

### Step 7.2 — Launch EC2 instance for LexicMap search

**Script:** `scripts/run_atb_lexicmap_ec2.sh`

```bash
bash scripts/run_atb_lexicmap_ec2.sh
```

This script:
1. Creates an SSH key pair (`~/.ssh/atb_search_key`) and imports to EC2
2. Creates a security group (SSH only)
3. Verifies/auto-detects the latest Amazon Linux 2023 ARM64 AMI
4. Launches an `r7g.4xlarge` on-demand instance (16 vCPU, 128 GB RAM)
   — ARM64 to match the LexicMap binary and minimise cost (~$0.80/hr)
5. User-data installs: `mountpoint-s3`, `LexicMap v0.8.1 (ARM64)`
6. Mounts the ATB LexicMap index from S3:
   ```
   mount-s3 --read-only --prefix 202408/ allthebacteria-lexicmap /mnt/atb_index --no-sign-request
   ```
7. Waits for instance ready (polls for `/home/ec2-user/READY`)
8. SCPs `kl_failing_query.fasta` to the instance
9. Runs LexicMap search:
   ```bash
   lexicmap search \
       -d /mnt/atb_index \
       /home/ec2-user/query.fasta \
       -o /home/ec2-user/results.tsv.gz \
       --align-min-match-pident 90 \
       --min-qcov-per-genome 70 \
       --top-n-genomes 0 \
       -j 16
   ```
10. Waits for `DONE` flag, fetches results, terminates instance

**Check progress:**
```bash
ssh -i ~/.ssh/atb_search_key -o StrictHostKeyChecking=no ec2-user@<PUBLIC_IP> \
    'tail -20 /home/ec2-user/lexicmap_run.log'

ssh -i ~/.ssh/atb_search_key -o StrictHostKeyChecking=no ec2-user@<PUBLIC_IP> \
    'test -f /home/ec2-user/DONE && echo DONE || echo Still running'
```

**Fetch when complete:**
```bash
bash scripts/fetch_atb_results.sh <INSTANCE_ID> <PUBLIC_IP>
```

> **Troubleshooting notes recorded:**
> - Spot quota = 0 on new AWS accounts → switch to on-demand
> - Default vCPU limit = 16 → use `r7g.4xlarge` not `c7g.8xlarge`
> - User-data 25,600 byte limit → do not embed large FASTA; SCP after boot
> - FUSE mount owned by root → run LexicMap via `sudo nohup bash`
> - Log must be in `/home/ec2-user/` not `/var/log/` (permission denied)

### Step 7.3 — Analyse LexicMap results

```bash
python3 scripts/analyse_atb_results.py
```

**Result:** 16,324,589 alignment rows.
- KL300: 329,240 candidate assemblies; top 20 at 99.8% qcov / 100% pident
- KL303: 221,599 candidate assemblies; top 5 at 99.8% qcov / 100% pident

Top candidates selected:

| Locus | Assembly | qcov% | pident% |
|-------|----------|-------|---------|
| KL300 | SAMD00053151 | 99.8 | 100.0 |
| KL300 | SAMEA103908486 | 99.8 | 100.0 |
| KL300 | SAMN12532904 | 99.8 | 100.0 |
| KL303 | SAMEA6656333 | 99.8 | 100.0 |
| KL303 | SAMEA6656328 | 99.8 | 100.0 |
| KL303 | SAMN20959615 | 99.8 | 100.0 |

### Step 7.4 — Validate candidates on M3

**Script:** `scripts/validate_kl300_kl303_m3.sh`

```bash
scp scripts/validate_kl300_kl303_m3.sh \
    ebenezef@m3.massive.org.au:~/EC-K-typing-G1G4/scripts/

# On M3:
sbatch ~/EC-K-typing-G1G4/scripts/validate_kl300_kl303_m3.sh
```

Extracts the 6 candidate FASTAs from ATB archives on M3, then runs
Kaptive normalised scoring. Output:
`DB/kl300_303_validation_summary.tsv`

**Monitor:**
```bash
tail -f ~/EC-K-typing-G1G4/logs/val_kl300_303_JOBID.log
```

### Step 7.5 — Fetch validated FASTAs and build v0.9 locally

```bash
# Fetch validation summary and FASTAs from M3
scp ebenezef@m3.massive.org.au:~/EC-K-typing-G1G4/DB/kl300_303_validation_summary.tsv \
    DB/

mkdir -p DB/kl300_303_fastas
scp "ebenezef@m3.massive.org.au:~/EC-K-typing-G1G4/DB/kl300_303_fastas/*.fa" \
    DB/kl300_303_fastas/
```

### Step 7.6 — Extract loci and build v0.9

**Script:** `scripts/replace_kl300_kl303.py`

```bash
# Auto-selects first correctly self-typed assembly per locus:
python3 scripts/replace_kl300_kl303.py --auto

# Or specify manually:
python3 scripts/replace_kl300_kl303.py \
    --kl300_fasta DB/kl300_303_fastas/SAMD00053151.fa \
    --kl303_fasta DB/kl300_303_fastas/SAMEA6656333.fa
```

The script:
1. Loads existing KL300/KL303 sequences from v0.8 as alignment queries
2. Uses minimap2 (`asm5` preset) to locate each locus on the validated assembly
3. Extracts the locus region ±500 bp
4. Re-annotates with pyrodigal + BLASTp
5. Replaces KL300/KL303 records in v0.8 → writes `DB/EC-K-typing_group1and4_v0.9.gbk`

### Step 7.7 — Re-run full validation on v0.9

```bash
# Update DB path in type_normalized.py to point to v0.9, then:
python3 scripts/type_normalized.py

# Expected: 93/93 correctly self-typed, 100% typeability
```

---

## Phase 8 — v0.9 Cleanup (Post-KL300/KL303 replacement)

After replacing KL300/KL303 and running full validation (651 loci × 1,126 assemblies),
the following additional issues were identified and resolved:

### Step 8.1 — Remove oversized novel loci (>60 CDS)

**Script:** `scripts/remove_oversized_loci.py`

12 novel ATB loci had >60 CDS because the extraction pipeline captured extensive
flanking chromosomal sequence. These universally scored ~2.0 against every assembly
(shared conserved chromosomal genes present in all *E. coli*), overriding correct
type assignments.

```bash
python scripts/remove_oversized_loci.py --dry-run   # preview
python scripts/remove_oversized_loci.py              # apply
```

Threshold: `--cds-threshold 60` (default). The original 93 BSI loci have a maximum
of 53 CDS (KL308), cleanly separating legitimate loci from extraction artefacts.

### Step 8.2 — Remove KL337

KL337 (BSI locus) was removed because:
- The cluster has only 2 members
- Both members self-type to other loci (KL389 and KL767 respectively)
- No assembly in the 1,126-genome validation set types uniquely to KL337
- There is no discriminating public representative

Removed from `DB/EC-K-typing_group1and4_v0.9.gbk`, `DB/KL_G1G4_mapping.tsv`,
and `DB/KL_G1G4_mapping_filtered.tsv`.

### Step 8.3 — Merge indistinguishable locus pairs

Three pairs of loci were identified where all genes of both loci were found at
100% identity in each other's representative assemblies (AS_norm = 2.0 for both
when scored against either assembly):

| Pair | Resolution | Rationale |
|------|------------|-----------|
| KL443 / KL446 | Merged → KL443 | KL443's rep self-types correctly; KL446's does not |
| KL716 / KL968 | Merged → KL716 | KL716's rep self-types correctly; KL968's does not |
| KL830 / KL843 | Merged → KL830 | KL830 has more genes (20 vs 19); wins AS tiebreaker |

Implementation: removed the subordinate locus from the GenBank, then remapped
its representative accession to the dominant locus in `novel_rep_kl_map.tsv`
and `novel_kl_summary.tsv`.

### Step 8.4 — Fix tiebreaker (see Normalised Scoring section)

Replaced 0.5% window + genes_found tiebreaker with pure AS_norm descending,
then raw AS descending. Updated in `scripts/type_normalized.py`.

### Step 8.5 — Final validation

```bash
python scripts/type_normalized.py \
    --db DB/EC-K-typing_group1and4_v0.9.gbk \
    --suffix v0.9norm_full \
    --threads 8
```

**Result: 651/651 (100%) self-typing; 1,126/1,126 (100%) typeability**

---

## Kaptive GenBank Format Requirements

Each record in the database must follow Kaptive's expected format:

```
LOCUS       KL392_ATB    45123 bp    DNA
...
FEATURES
  source    1..45123
            /note="K locus: KL392"
            /organism="Escherichia coli"
            /mol_type="genomic DNA"
  CDS       <start>..<end>
            /gene="wzy"
            /product="Wzy polymerase"
            /translation="MAST..."
```

The `source` feature **must** contain `/note="K locus: KLxxx"` — this is
how Kaptive identifies the locus name. The `name` field of the SeqRecord
is also set to `KLxxx_ATB` for uniqueness.

---

## Normalised Scoring — Technical Details

**Script:** `scripts/normalise_novel_scores.py`

Raw Kaptive `--scores` output provides, per assembly per locus:
- `AS`: total alignment score
- `genes_found`: number of expected CDS with alignments
- `genes_expected`: total expected CDS in that locus
- `q_len`: query (locus) length

Normalised score:
```
AS_norm = AS / total_expected_gene_bp
```

where `total_expected_gene_bp` is the sum of CDS lengths for that locus
(computed from the GenBank). If unavailable, `q_len` is used as a fallback.

Locus assignment (v0.9 tiebreaker):
1. Compute `AS_norm` for all loci with `AS > 0`
2. Sort descending by `AS_norm`, then descending by raw `AS` as tiebreaker
3. Select the top-ranked locus
4. Mark as **Typeable** if `genes_found / genes_expected >= 0.50`

> **Note (v0.9 change):** Prior versions used a 0.5% window + `genes_found` tiebreaker
> (steps 2–4: shortlist within `AS_norm >= top * 0.995`, then pick highest `genes_found`).
> This was replaced with pure `AS_norm` descending then raw `AS` descending, which
> eliminates systematic bias toward larger loci in near-tie cases and resolved the
> KL443/KL446, KL716/KL968, and KL830/KL843 indistinguishable-pair failures.

---

## Version History

| Version | Loci | Notes |
|---------|------|-------|
| v0.1–v0.3 | ~45 | Initial BSI isolate set, various fixes |
| v0.4 | 63 | Representative updates |
| v0.5 | 78 | Additional loci added |
| v0.6 | 93 | Final BSI-only set; KL300–KL391 + K24 + K96 |
| v0.7 | — | Deprecated: combined all-groups DB (wzy-interference artefact) |
| v0.8 | 667 | + 574 ATB-derived novel loci (KL392–KL968); KL486/562/812 removed |
| **v0.9** | **651** | **KL300/KL303 replaced; 15 oversized/artefact loci removed; KL337 removed; 3 indistinguishable pairs merged; tiebreaker fixed; 651/651 (100%) self-typing** |

---

## Key Design Decisions

1. **G1/G4 only, not combined:** Combined all-groups DB deprecated at v0.7
   due to wzy-interference. Sequential typing (G2/G3 first) is mandatory.

2. **KL numbering:** New loci start at KL392 (not KL344 as the clustering
   script comments suggest) because the v0.6 database already contains loci
   up to KL391.

3. **Normalised scoring over raw AS:** Raw Kaptive scores are proportional
   to locus size. Normalisation by total expected CDS length is essential
   for correct typing when loci of very different sizes share gene content
   (e.g., partial captures).

4. **Bidirectional clustering (--cov-mode 0):** Ensures that both the query
   and the target must meet the 80% coverage threshold. This prevents short
   sequences from clustering into larger, unrelated loci. Note: partial
   captures shorter than the reference can still slip through if the
   reference has low coverage of the novel sequence — handled by the
   additional BLASTn partial-capture screening step in the annotation script.

5. **Representative selection:** Cluster representatives in MMseqs2 are
   chosen by the algorithm (typically the longest or most connected sequence).
   For KL300/KL303, the ATB representatives were replaced using LexicMap
   to find assemblies with near-perfect locus matches.

6. **Minimum locus length:** Novel extractions < 15,000 bp are discarded
   as likely partial captures due to contig breaks within the locus.
