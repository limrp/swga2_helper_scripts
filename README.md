# SWGA2 primer analysis helper scripts

Helper scripts to analyze **primer sets selected by [swga2](https://github.com/songlab-cal/swga2)** for Selective Whole Genome Amplification (SWGA).

These tools are designed to:

- Inspect **genome-wide spacing** of primer binding sites (both strands).
- Quantify **gap statistics** (mean/max distances between sites).
- Specifically characterize **converging forward→reverse** (F→R) pairs, which are most relevant for Φ29-based SWGA.
- Produce **plots** and **summary tables** to compare primer sets and decide which ones to use.

The repository is intended to accompany swga2’s output (e.g. `step4_df.csv`) and give an extra level of **quality control** and **biological interpretability** for primer sets, both on plasmids and on full genomes (including multi-chromosome genomes such as *Toxoplasma gondii*).

---

## Contents

Currently the main helper script is:

- `summarize_gaps.py`  
  Analyze one or more primer sets against a reference genome (plasmid or genome, linear or circular), compute spacing metrics, and generate plots.

More scripts can be added later (for example: GC bias inspection, off-target summaries, or integration with swga2’s HDF5 metrics).

---

## Installation

### Requirements

- Python 3.8+
- Recommended packages:
  - `pandas`
  - `matplotlib`

You can install them with:

```bash
pip install pandas matplotlib
```

or create a small virtual environment:

```bash
python -m venv venv
source venv/bin/activate      # Linux/macOS
# .\venv\Scripts\activate     # Windows PowerShell

pip install pandas matplotlib
```

---

## Input formats

### 1. Primer sets from swga2 (`step4_df.csv` style)

The typical swga2 output for top primer sets (`step4_df.csv`) has lines like:

```text
[AATTAATTC, ACGGTAAAT, ACGTCAATA, ATAAACCAG, CAACTACAA, CACTCAAAG, CATATATGG], 49130.9
[AATAAACCA, AATTAATTC, ACGGTAAAT, ACGTCAATA, CAACTACAA, CACTCAAAG, CATATATGG], 49130.9
...
```

- Each line is:
  - a primer set in brackets: `[P1, P2, ..., Pk]`
  - followed by a comma and a **score** (e.g. predicted on-target amplification, as defined by swga2).

Use `--step4-format` to parse this style correctly.

### 2. Generic CSV format

Alternatively, a generic CSV is supported:

```text
AATTAATTC,ACGGTAAAT,ACGTCAATA,ATAAACCAG,CAACTACAA,CACTCAAAG,CATATATGG,49130.9
AATAAACCA,AATTAATTC,ACGGTAAAT,ACGTCAATA,CAACTACAA,CACTCAAAG,CATATATGG,49130.9
...
```

- First column (or row[0]) is interpreted as **comma-separated primers**.
- Optional second column is interpreted as a numeric **score**.
- Use this mode **without** `--step4-format`.

### 3. Genome FASTA

You must provide a reference sequence in FASTA:

- Can be a **single plasmid** (circular or linear).
- Or a **multi-FASTA genome** (e.g., multiple chromosomes for *T. gondii*).

Headers are used only for naming; sequences are concatenated per record and converted to uppercase.

Example:

```text
>pcDNA
AGCTGATCGATCGATCGATC...
```

or

```text
>chr1_Tgondii
...
>chr2_Tgondii
...
...
```

---

## `summarize_gaps.py`

### Purpose

Given:

- A set of SWGA primers (one or many sets).
- A reference genome (plasmid or multi-chromosome genome).
- A flag indicating whether each sequence should be treated as **circular** or **linear**.

`summarize_gaps.py` will:

1. Find **all binding sites** for each primer and its reverse complement in each sequence.
2. Compute **gaps between consecutive sites** along each chromosome/plasmid:
   - **All gaps** (ignoring strand).
   - **Converging forward→reverse gaps** (F→R only), i.e., plus-strand site followed by downstream minus-strand site.
3. Pool gaps across chromosomes.
4. Output:
   - A **CSV summary** of metrics *per primer set*.
   - **Plots** per primer set: histograms and sorted-gap plots.

### Command-line usage

```bash
python summarize_gaps.py \
  -f step4_df.csv \
  -g genome.fasta \
  -o output_dir \
  -c true \
  --step4-format
```

#### Arguments

- `-f`, `--input_file`  
  Path to the primer file.  
  - Use `--step4-format` if it’s swga2-style (`[P1, ..., Pk], score`).
  - Otherwise, a generic CSV with comma-separated primers (and optional score).

- `-g`, `--input_genome`  
  Genome FASTA (plasmid or multi-FASTA).

- `-o`, `--output_dir`  
  Output directory where CSV and plots will be written.  
  - The parent directory **must exist**.
  - The script will create the final directory if needed.

- `-c`, `--circular` (`true` or `false`)  
  Treat each sequence as circular (e.g. plasmids) or linear (e.g. chromosomes).

- `--step4-format`  
  Indicate that the input file is swga2’s `step4_df.csv` format:
  - Lines like `[P1, P2, ..., Pk], score`
  - No header line.

> **Note**: if you add later a parameter for template length (e.g. `--max_amplicon_length`), you can use it to filter F→R gaps that are physically plausible on single molecules.

---

## Output

### 1. Summary metrics (`gap_stats_per_set.csv`)

The script writes a CSV file in `output_dir`:

```text
gap_stats_per_set.csv
```

with one row per primer set, e.g.:

```csv
set_index,score,n_primers,n_chromosomes,n_chrom_with_sites,n_sites,n_gaps_all,mean_gap_all,max_gap_all,n_gaps_fr,mean_gap_fr,max_gap_fr
1,49130.9,7,13,13,845,832,73250.4,210345,412,62890.3,185210
2,49130.9,7,13,13,823,810,75510.8,198765,401,64532.7,176432
...
```

**Columns**

- `set_index`  
  1-based index of the primer set in the input file.

- `score`  
  The swga2 score for this set (or NaN if not provided).

- `n_primers`  
  Number of primers in this set.

- `n_chromosomes`  
  Number of sequences in the FASTA.

- `n_chrom_with_sites`  
  How many sequences have at least one binding site for this set.

- `n_sites`  
  Total number of binding sites (all chromosomes, both strands).

- `n_gaps_all`  
  Count of all gaps between consecutive binding sites (ignoring strand).

- `mean_gap_all`, `max_gap_all`  
  Mean and max distance (bp) between any consecutive binding sites (ignoring direction).

- `n_gaps_fr`  
  Number of **forward→reverse** converging gaps (plus strand site followed by downstream minus strand).

- `mean_gap_fr`, `max_gap_fr`  
  Mean and max distance (bp) between converging F→R pairs.

> **Biological interpretation (SWGA context)**  
> - **Gaps_all** describe generic primer density and uniformity over the genome.  
> - **F→R gaps** approximate “amplicon-like” spans where a forward primer and a downstream reverse primer can cooperate to strongly amplify a region with Φ29.  
> - A small `max_gap_fr` (relative to your acceptable SWGA amplicon scale) suggests no large “deserts” lacking converging primer pairs.

---

## Acknowledgements

- Primer set generation: [swga2](https://github.com/songlab-cal/swga2) (Song Lab).  
- This repository focuses on **post-selection analysis** of those primer sets: spacing, F→R convergence, and plotting.



