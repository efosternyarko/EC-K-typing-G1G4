# EC-K-typing Group 1 & 4

A [Kaptive](https://github.com/klebgenomics/Kaptive)-compatible reference database for typing *Escherichia coli* **Group 1 and Group 4** capsular K-loci.

**Current version: v1.2** — 623 reference loci · 100% self-typing · 100% typeability

This database complements [EC-K-typing](https://github.com/rgladstone/EC-K-typing) (Group 2 & 3, Gladstone et al.) to enable K-locus typing across all four *E. coli* capsule groups.

---
## Prerequisites

- [Kaptive](https://github.com/klebgenomics/Kaptive) ≥ 3.0 (`pip install kaptive`)
- Python ≥ 3.8 with `pandas` and `biopython` (`pip install pandas biopython`)

---

## Background

*E. coli* capsule loci fall into four biosynthetic groups:

| Group | Pathway | Database |
|-------|---------|----------|
| 1 | Wzy-dependent (*cps* locus) | **This database** |
| 2 | ABC transporter (*kps* locus) | [EC-K-typing (Gladstone)](https://github.com/rgladstone/EC-K-typing) |
| 3 | ABC transporter (*kps* locus) | [EC-K-typing (Gladstone)](https://github.com/rgladstone/EC-K-typing) |
| 4 | Wzy-dependent (*cps* locus) | **This database** |

**G1/G4 and G2/G3 are mutually exclusive** — no *E. coli* isolate carries both. They must be typed in two separate steps, not together. Running both databases simultaneously causes a *wzy*-interference artefact that misclassifies genuine G2/G3 isolates as untypeable.

---

## Workflow

### Optional pre-screen: kpsM marker (faster alternative to Step 1)

For large collections where G2/G3 Kaptive runtime is prohibitive, you can
pre-screen assemblies for the `kpsM` gene — an inner-membrane ABC-transporter
component present in all Group 2/3 strains and absent from Group 1/4 strains.
Assemblies that are kpsM-negative can be passed directly to Step 2.

> **Note:** The kpsM screen identifies Group 2/3 organisms but does **not**
> assign a K-locus type. Use it only when you do not need G2/G3 K-locus
> assignments. For full G2/G3 typing, use Step 1 as described below.

**Requirements:** `minimap2` on your `PATH` (`conda install -c bioconda minimap2` or `pip install minimap2`).

```bash
# Screen a directory of assemblies; run Kaptive on kpsM-negative ones:
python3 scripts/kpsm_screen.py \
    --assembly-dir assemblies/ \
    --kpsm-ref DB/kpsM_reference.fasta \
    --kaptive-db DB/EC-K-typing_group1and4_v1.2.gbk \
    --output-dir results/kpsm_screen/ \
    --run-kaptive

# Screen only (skip Kaptive, inspect kpsm_screen_results.tsv first):
python3 scripts/kpsm_screen.py \
    --assembly-dir assemblies/ \
    --kpsm-ref DB/kpsM_reference.fasta \
    --output-dir results/kpsm_screen/
```

**Output:** `results/kpsm_screen/kpsm_screen_results.tsv` — one row per
assembly with columns: `assembly`, `kpsM_hit` (True/False), `kpsM_pident`,
`kpsM_qcov`, `group` (`Group2_3` or `G1_G4_candidate`).

kpsM-positive assemblies (`group = Group2_3`) are Group 2/3 strains — G1/G4
typing is not applicable to them. kpsM-negative assemblies
(`group = G1_G4_candidate`) should be passed to Step 2 (or Step 3 directly
if you have used `--run-kaptive`).

If you use this pre-screen, skip Step 1 and proceed to Step 2 using only
the kpsM-negative assemblies.

---

### Step 1 — G2/G3 typing

Run the Gladstone Group 2 & 3 database on all your assemblies (replace `assemblies/` with the path to your assembly folder and database):


```bash
mkdir -p results_G23
kaptive assembly \
     DB/EC-K-typing_group2and3_v3.0.0.gbk \
     assemblies/*.fasta \
    -o results_G23/kaptive_results.tsv
```
Here Kaptive is performing typing of assemblies with `assembly` parameter. All other parameters are set to the default.

**Example output** (`kaptive_results.tsv`, key columns shown):

```
Assembly    Best match locus  Best match type  Match confidence  Problems  Identity  Coverage
sample_A    KL2               KL2              Perfect                     100.0%    100.0%
sample_B    KL127             KL127            Very High         !         99.8%     98.4%
sample_C    none              none             Untypeable
sample_D    KL15              KL15             Perfect                     100.0%    100.0%
sample_E    none              none             Untypeable
```

Assemblies with `Match confidence = Untypeable` have no G2/G3 locus — these are your G1/G4 candidates and should be passed to Step 2.

---

### Step 2 — G1/G4 typing

#### Collect untypeable assemblies from Step 1

```bash
# Extract the names of untypeable assemblies
grep -w "Untypeable" results_G23/kaptive_results.tsv | cut -f1 > untypeable_ids.txt

# Copy their FASTA files to a new directory
mkdir -p untypeable/
while IFS= read -r id; do
    cp "assemblies/${id}.fasta" untypeable/
done < untypeable_ids.txt
```

#### Run Kaptive in scores mode on the untypeables


```bash
mkdir -p results_G14
kaptive assembly \
    DB/EC-K-typing_group1and4_v1.2.gbk \
    untypeable/*.fasta \
    --scores results_G14/kaptive_scores.tsv \
    -t 8
```

The `--scores` flag produces a full locus × assembly score matrix. This intermediate file is required for the normalisation step — it is **not** the final result.

**Example scores matrix** (`results_G14/kaptive_scores.tsv`, first few rows):

```
Assembly    Locus   AS      mlen    blen    q_len   genes_found  genes_expected
sample_C    KL300   45231   45380   38904   45380   36           36
sample_C    KL302   68207   54210   45992   54210   38           40
sample_C    KL304   12445   52546   11230   52546   10           44
sample_C    KL436   58912   30210   28740   30210   16           16
...         ...     ...     ...     ...     ...     ...          ...
sample_E    KL300   22104   45380   19023   45380   18           36
sample_E    KL302   31882   54210   27410   54210   20           40
...         ...     ...     ...     ...     ...     ...          ...
```

> This file contains one row per assembly per reference locus (623 loci × number of assemblies). The `AS` column is the raw alignment score — **not** directly interpretable without normalisation. Step 3 converts this into a per-locus best-match assignment with a meaningful score.

---

### Step 3 — Normalise scores

The raw `AS` values in the scores matrix are not comparable across loci of different sizes — a locus with 50 genes will always score higher than one with 10 genes, even at the same per-base identity. Normalisation divides AS by the total expected CDS length of each reference locus, producing a size-independent metric (Norm AS) where ~2.0 = perfect match.

Download [`scripts/normalise_kaptive_scores.py`](scripts/normalise_kaptive_scores.py) from this repository (see [Database files](#database-files) below), then run it to convert raw scores to normalised per-base identity values (see run command below).

#### Run command

```bash
python3 scripts/normalise_kaptive_scores.py \
    --db DB/EC-K-typing_group1and4_v1.2.gbk \
    --in results_G14/kaptive_scores.tsv \
    --out results_G14/kaptive_results_norm.tsv
```

> `python3` is recommended. If Python 3 is already your system default, `python` works equally well for the command above.

**Example output** (`kaptive_results_norm.tsv`):

```
Assembly    Best match locus  Best match confidence  Genes found  Genes expected  Gene coverage  Raw AS  Norm AS
sample_C    KL302             Typeable               38           40              95.0%          68207   1.8059
sample_D    KL304             Typeable               44           44              100.0%         93072   2.0
sample_E    KL306             Typeable               37           37              100.0%         75810   2.0
sample_F    KL436             Typeable               16           16              100.0%         29576   1.9545
sample_G    KL305             Typeable               50           50              100.0%         88182   2.0
sample_H    none              Untypeable             3            22              13.6%          5204    0.2741
```

---

## Interpreting results

### Column guide

| Column | What it means |
|--------|---------------|
| `Assembly` | Filename stem as passed to Kaptive |
| `Best match locus` | Assigned K-locus (e.g. `KL302`, `K24`). `none` = untypeable |
| `Best match confidence` | `Typeable` if ≥ 50% of expected genes found; `Untypeable` otherwise |
| `Genes found` / `Genes expected` | Number of reference CDS detected vs total expected |
| `Gene coverage` | `Genes found ÷ Genes expected` as a percentage |
| `Raw AS` | Cumulative alignment score — **not** comparable across loci of different sizes (range: ~10,000–100,000+) |
| `Norm AS` | `Raw AS ÷ total expected CDS length` — size-independent, comparable across all loci (~2.0 = perfect) |

### Interpreting Norm AS

| Norm AS | Meaning |
|---------|---------|
| ~2.00 | Perfect match: 100% identity across all reference genes |
| ≥ 1.90 | High confidence: all or nearly all genes found at high identity |
| 1.50 – 1.90 | Moderate: partial locus or some sequence divergence |
| < 1.50 | Low confidence: possibly novel, highly divergent, or fragmented |

### Why normalised scoring?

Standard Kaptive locates a locus by first aligning the full reference sequence to the assembly. If this step fails — because the locus is fragmented across contigs or too divergent — the assembly is marked **Untypeable**, even though individual locus genes are present. This affects ~57% of G1/G4 assemblies.

The `--scores` + normalisation approach instead aligns each reference gene individually against the whole assembly, succeeding even on fragmented or divergent assemblies, achieving **100% typeability**.

| Mode | Typeability | Notes |
|------|-------------|-------|
| Standard Kaptive | ~43% for G1/G4 | Also gives locus position and structure information |
| Scores + normalisation | **100%** | Recommended; use for final type assignments |

---

## Database files

Only these files are needed for typing:

| File | Description |
|------|-------------|
| `DB/EC-K-typing_group1and4_v1.2.gbk` | **G1/G4 database — 623 loci** |
| `DB/EC-K-typing_group2and3_v3.0.0.gbk` | Gladstone G2/G3 database (included for convenience) |
| `DB/kpsM_reference.fasta` | kpsM/kpsM_3 reference sequences for optional Group 2/3 pre-screen |
| `scripts/normalise_kaptive_scores.py` | Normalisation script for Step 3 |

---

## Locus nomenclature

| Name | Description |
|------|-------------|
| `K24` | Classical K-type designation |
| `KL300`–`KL391` | Novel loci from *E. coli* BSI isolates (EnteroBase) |
| `KL392`–`KL967` | Novel loci from AllTheBacteria (~2.4 M assemblies) |

Most G1/G4 loci do not yet have a classical K-antigen serotype name. The KL designation is both the locus name and the type assignment reported by Kaptive.

---

## Citation

If you use this database, please cite:

- This repository (https://github.com/efosternyarko/EC-K-typing-G1G4)
- Gladstone RA et al. (2026). *E. coli* capsule typing. [*Nature Microbiology*](https://doi.org/10.1038/s41564-026-02283-w)
- [Kaptive](https://github.com/klebgenomics/Kaptive): Lam MMC et al. (2022). [doi:10.1099/mgen.0.000800](https://doi.org/10.1099/mgen.0.000800)

## Contributors

- [Emmanuel Foster-Nyarko](https://github.com/efosternyarko)
- [Louise Cerdeira](https://github.com/lcerdeira)
- Mary Maranga

## License

MIT License — see [LICENSE](LICENSE).
