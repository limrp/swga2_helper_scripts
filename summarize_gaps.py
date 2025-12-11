#!/usr/bin/env python3
"""
Summarize spacing between primer binding sites for multiple primer sets
over a (possibly multi-chromosome) genome.

For each primer set:
- Find binding sites (both strands) in each chromosome.
- Compute gaps between consecutive binding sites (per chromosome).
- Pool all gaps across chromosomes to compute global statistics.
- Generate global plots (all gaps, all vs F/R gaps, sorted gaps).

Usage example (step4_df style input):

    python summarize_gaps.py \
        -f step4_df.csv \
        -g Tgondii_TGA4.fasta \
        -o tgondii_gaps \
        -c false \
        --step4-format
"""

import argparse  # For command-line interface
import math  # For NaN, mean, etc.
from pathlib import Path  # For filesystem paths in an OS-independent way

import matplotlib.pyplot as plt  # For plots
import pandas as pd  # For generic CSV reading

# ===========================
# CLI helpers
# ===========================


def str_to_bool(v):
    """
    Convert a command-line string to a boolean.

    Allows values like: "yes", "no", "true", "false", "1", "0".
    Raises an error if the string cannot be interpreted as boolean.
    """
    # If it's already a boolean, just return it.
    if isinstance(v, bool):
        return v

    # Convert to lowercase for easier comparison.
    val = v.lower()

    # Map common "true" variants to True.
    if val in ("yes", "y", "true", "t", "1"):
        return True

    # Map common "false" variants to False.
    if val in ("no", "n", "false", "f", "0"):
        return False

    # If not recognized, raise an argparse error.
    raise argparse.ArgumentTypeError(f"Invalid boolean value: '{v}'")


def parse_args():
    """
    Define and parse the command-line arguments.

    Returns an argparse.Namespace object with all parsed options.
    """
    # Create the top-level parser with a description.
    parser = argparse.ArgumentParser(
        description="Summarize primer sets with spacing metrics over a genome."
    )

    # Input primer file (step4_df style or generic CSV).
    parser.add_argument(
        "-f",
        "--input_file",
        required=True,
        type=Path,  # Convert the string path directly into a Path object
        help=(
            "Input primer file. If --step4-format is set, expect step4_df.csv style "
            "lines: '[p1, p2, ...], score' (no header). Otherwise, expect a CSV where "
            "the first column contains comma-separated primers and an optional second "
            "column contains a numeric score."
        ),
    )

    # Genome FASTA (can be multi-FASTA with multiple chromosomes).
    parser.add_argument(
        "-g",
        "--input_genome",
        required=True,
        type=Path,
        help="Genome FASTA. Can be multi-FASTA (e.g. multiple chromosomes).",
    )

    # Output directory for summary CSV + plots.
    parser.add_argument(
        "-o",
        "--output_dir",
        required=True,
        type=Path,
        help="Output directory for summary CSV and plots.",
    )

    # Flag for circular vs linear sequences (applied to each sequence).
    parser.add_argument(
        "-c",
        "--circular",
        type=str_to_bool,
        default=False,
        help="Treat each sequence (plasmid/chromosome) as circular. Default False.",
    )

    # Flag to indicate step4_df-like format for primers.
    parser.add_argument(
        "--step4-format",
        action="store_true",  # Becomes True when present
        help="Input is step4_df.csv style: '[p1, p2, ...], score' with no header.",
    )

    # Max F→R distance to consider in the SWGA2-style F→R metric.
    parser.add_argument(
        "--within-threshold",
        type=int,
        default=100000,
        help=(
            "Max forward→reverse gap to include (bp) for the SWGA2-style F→R metric. "
            "Matches SWGA2's get_positional_gap_lengths_alternating threshold. "
            "Default 100000."
        ),
    )

    # Parse the arguments from the command line and return them.
    return parser.parse_args()


# ===========================
# FASTA + primer parsing
# ===========================


def load_fasta_dict(path: Path):
    """
    Load a (multi-)FASTA file into a dictionary: {seq_name: sequence}.

    - 'seq_name' is taken from the header line (everything after '>')
      up to the first whitespace.
    - The sequence is returned as an uppercase string with newlines removed.
    """
    sequences = {}  # Will store {name: sequence}
    name = None  # Current sequence name
    seq_chunks = []  # List of lines for the current sequence

    # Open the FASTA file for reading.
    with path.open() as f:
        for line in f:
            line = line.strip()  # Remove trailing newline and spaces

            # Skip empty lines (if any).
            if not line:
                continue

            if line.startswith(">"):
                # This is a header line, so we start a new sequence.
                # First, if we already had a previous sequence, save it.
                if name is not None:
                    sequences[name] = "".join(seq_chunks).upper()

                # Now parse the new name (without the '>' and trimming at first space).
                name = line[1:].split()[0]
                seq_chunks = []  # Reset chunks for the new sequence
            else:
                # This is a sequence line; accumulate it.
                seq_chunks.append(line)

    # After the loop ends, we might have one last sequence to store.
    if name is not None:
        sequences[name] = "".join(seq_chunks).upper()

    # Return a dict mapping sequence names to their sequences.
    return sequences


def parse_step4_line(line: str):
    """
    Parse a step4_df-style line of the form:

        [AATTAATTC, ACGGTAAAT, ... , CATATATGG], 49130.9

    Returns:
       primers (list of strings), score (float)
    """
    # Remove leading/trailing whitespace.
    line = line.strip()

    # If the line is empty, we return (None, None) as a signal.
    if not line:
        return None, None

    # We split at the last occurrence of '],' to separate primers and score.
    try:
        primers_part, score_part = line.rsplit("],", 1)
    except ValueError:
        # If rsplit fails, the format is not as expected.
        raise ValueError(f"Cannot parse line (step4-format expected): {line}")

    # Clean trailing/leading spaces from the part containing primers.
    primers_part = primers_part.strip()

    # Remove leading '[' if present, leaving "P1, P2, P3, ..."
    if primers_part.startswith("["):
        primers_part = primers_part[1:]

    # Split by comma and strip each primer.
    primers = [p.strip() for p in primers_part.split(",") if p.strip()]

    # Remaining part should be the score string, e.g. "49130.9".
    score = float(score_part.strip())

    return primers, score


def load_primer_sets_step4(path: Path):
    """
    Load primer sets from a step4_df.csv-style file.

    Each non-empty line should be of the form:
        [P1, P2, ..., Pk], score

    Returns:
        list of (primers_list, score_float)
    """
    primer_sets = []  # Will collect (primers, score) for each line

    # Open the file and iterate line by line.
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                # Skip empty lines.
                continue
            primers, score = parse_step4_line(line)
            primer_sets.append((primers, score))

    return primer_sets


def load_primer_sets_generic_csv(path: Path):
    """
    Load primer sets from a generic CSV file (no header).

    Assumptions:
      - First column: comma-separated primers (e.g. 'P1,P2,P3').
      - Optional second column: numeric score for that primer set.

    Returns:
      list of (primers_list, score_or_nan)
    """
    # Read the CSV with no header, so columns are indexed as 0, 1, ...
    df = pd.read_csv(path, header=None)

    primer_sets = []

    # Iterate over each row in the DataFrame.
    for _, row in df.iterrows():
        # Convert column 0 to string, then split by comma and strip each primer.
        primers_str = str(row[0])
        primers = [p.strip() for p in primers_str.split(",") if p.strip()]

        # If a second column exists and is not NaN, use it as score, else NaN.
        score = float(row[1]) if len(row) > 1 and not pd.isna(row[1]) else math.nan

        primer_sets.append((primers, score))

    return primer_sets


# ===========================
# Binding sites + gap computation
# ===========================

# Translation table for reverse complement (A <-> T, C <-> G).
comp_table = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    """
    Return the reverse complement of a DNA sequence string s.
    """
    # translate() replaces each base by its complement; [::-1] reverses the string.
    return s.translate(comp_table)[::-1]


def find_sites(genome_seq: str, primer: str):
    """
    Find all binding sites for a primer and its reverse complement on one sequence.

    Returns a list of dicts of the form:
        {"pos": start_index, "strand": "+" or "-"}

    - Positions are 0-based indices (Python style).
    - "+" means the primer itself matches the forward strand.
    - "-" means the reverse complement matches (i.e. binding on the opposite strand).
    """
    # Normalize primer to uppercase.
    primer = primer.upper()

    # Compute reverse complement of the primer.
    rc = revcomp(primer)

    sites = []  # Will hold all found sites (forward + reverse)

    # ----- Forward strand matches -----
    start = 0  # Start scanning at position 0
    while True:
        # Find the next occurrence of primer starting at 'start'.
        pos = genome_seq.find(primer, start)
        if pos == -1:
            # .find returns -1 when no more matches.
            break
        # Record site with strand "+"
        sites.append({"pos": pos, "strand": "+"})
        # Move start one base ahead to allow overlapping matches if they exist.
        start = pos + 1

    # ----- Reverse strand matches -----
    start = 0  # Reset start for scanning the reverse complement
    while True:
        # Find the next occurrence of reverse complement rc.
        pos = genome_seq.find(rc, start)
        if pos == -1:
            break
        # Record site with strand "-"
        sites.append({"pos": pos, "strand": "-"})
        start = pos + 1

    return sites


def compute_gaps_from_positions(positions, L: int, circular: bool):
    """
    Given a list of binding site positions on one sequence, compute gaps
    between consecutive positions.

    - positions: list of integer indices (0-based).
    - L: length of the sequence.
    - circular: if True, include the wrap-around gap from last to first.

    Returns:
        list of gap lengths (integers, in bp).
    """
    # If we have fewer than 2 positions, we cannot form any gap.
    if len(positions) < 2:
        return []

    # Sort positions in ascending order along the sequence.
    positions = sorted(positions)

    gaps = []  # Will store gap lengths

    # Gaps between consecutive positions: pos[i+1] - pos[i]
    for i in range(len(positions) - 1):
        gaps.append(positions[i + 1] - positions[i])

    # If the sequence is circular, we also add the last->first gap:
    # L - last_pos + first_pos
    if circular:
        gaps.append(L - positions[-1] + positions[0])

    return gaps


def compute_fr_gaps(sites_sorted, L: int, circular: bool):
    """
    Compute gaps between *converging* forward→reverse site pairs.

    sites_sorted: list of dicts {"pos": int, "strand": "+" or "-"},
                sorted by genomic position.
    L: sequence length.
    circular: if True, include wrap-around pair (last → first) using modulo.

    Only keep gaps where:
        upstream site is "+" (forward) and downstream site is "-" (reverse),
    i.e. a converging forward→reverse pair.
    """
    n = len(sites_sorted)

    # If fewer than 2 sites, no gaps.
    if n < 2:
        return []

    gaps = []  # List to store F/R or R/F gap lengths

    # For circular sequences, inspect (i, (i+1) % n) for all i.
    # For linear sequences, inspect (i, i+1) for i = 0..n-2.
    if circular:
        indices = range(n)
    else:
        indices = range(n - 1)

    for i in indices:
        # Determine the index of the "next" site (j).
        j = (i + 1) % n if circular else (i + 1)

        # Get the two sites.
        s_i = sites_sorted[i]
        s_j = sites_sorted[j]

        # We only care about converging forward→reverse:
        # upstream "+" (forward), downstream "-" (reverse)
        if not (s_i["strand"] == "+" and s_j["strand"] == "-"):
            continue

        # Compute genomic distance from upstream F to downstream R
        if circular:
            # Use modulo to handle wrap-around.
            gap = (s_j["pos"] - s_i["pos"]) % L
        else:
            # For linear case, j > i so just subtract.
            gap = s_j["pos"] - s_i["pos"]

        gaps.append(gap)

        # DEBUG
        # DEBUG: now runs once for each F→R gap
        print(
            f"[DEBUG F→R] {s_i['pos']}({s_i['strand']}) -> {s_j['pos']}({s_j['strand']}), gap = {gap}"
        )

    return gaps


def compute_fr_gaps_alternating(positions_forward, positions_reverse, threshold: int):
    """
    SWGA2-style F→R gaps:
    For each forward site, find all reverse sites that are at the same or greater
    position and at most 'threshold' bp away, and record their distances.

    This mirrors the logic of optimize.get_positional_gap_lengths_alternating
    and the manual loop in summarize_sets.py.

    Parameters
    ----------
    positions_forward : list[int]
        Positions (0-based) of binding sites on the forward strand.
    positions_reverse : list[int]
        Positions (0-based) of binding sites on the reverse strand.
    threshold : int
        Maximum allowed distance (bp) between F and R to be included.

    Returns
    -------
    list[int]
        All forward→reverse distances (in bp) that satisfy:
        R >= F  and  R - F <= threshold.
    """
    # Ensure positions are sorted (SWGA2 assumes sorted lists).
    fwd = sorted(positions_forward)
    rev = sorted(positions_reverse)

    gaps = []
    i = j = 0

    # Two-pointer scan: walk forward and reverse lists together.
    while i < len(fwd) and j < len(rev):
        # Move the reverse pointer until it is not behind the current forward site.
        while j < len(rev) and rev[j] < fwd[i]:
            j += 1

        # Now try all reverse sites starting at 'j' that are within threshold.
        sub_j = j
        while (
            sub_j < len(rev)
            and rev[sub_j] - fwd[i] <= threshold
            and rev[sub_j] >= fwd[i]
        ):
            gaps.append(rev[sub_j] - fwd[i])
            sub_j += 1

        # Advance to the next forward site.
        i += 1

    return gaps


def safe_mean(lst):
    """
    Compute the mean of a list, or return NaN if the list is empty.
    """
    return sum(lst) / len(lst) if lst else math.nan


# def analyze_primer_set(primers, genome_dict, circular: bool):
#     """
#     Analyze one primer set over a (possibly multi-chromosome) genome.

#     - genome_dict: {chrom_name: chrom_seq}
#     - circular: treat each chromosome as circular if True, else linear.

#     For each chromosome:
#         - Find binding sites (forward + reverse).
#         - Compute gaps (all) and F/R gaps.

#     We do *not* connect chromosomes with gaps; each chromosome is independent.

#     Returns:
#         stats dict with keys:
#             n_primers, n_chromosomes, n_chrom_with_sites,
#             n_sites, n_gaps_all, mean_gap_all, max_gap_all,
#             n_gaps_fr, mean_gap_fr, max_gap_fr,
#             gaps_all, gaps_fr
#     """
#     # Total number of binding sites across all chromosomes.
#     total_sites = 0

#     # Total number of chromosomes in the genome.
#     n_chromosomes = len(genome_dict)

#     # Number of chromosomes that have at least one binding site.
#     n_chrom_with_sites = 0

#     # Global lists of gaps pooled across all chromosomes.
#     global_gaps_all = []
#     global_gaps_fr = []

#     # Iterate over chromosomes in the genome.
#     for chrom, seq in genome_dict.items():
#         # Length of this chromosome/sequence.
#         L = len(seq)

#         # Collect all binding sites for this primer set in this chromosome.
#         all_sites_chrom = []
#         for pr in primers:
#             # Find sites for each primer and accumulate them.
#             all_sites_chrom.extend(find_sites(seq, pr))

#         # If no sites on this chromosome, skip it.
#         if not all_sites_chrom:
#             continue

#         # We found sites in this chromosome.
#         n_chrom_with_sites += 1
#         total_sites += len(all_sites_chrom)

#         # ----- Compute "all gaps" on this chromosome -----
#         # Extract positions only, and sort them.
#         positions = sorted(s["pos"] for s in all_sites_chrom)

#         # Compute gaps between consecutive positions (and wrap if circular).
#         gaps_all_chrom = compute_gaps_from_positions(
#             positions,
#             L,
#             circular=circular,
#         )

#         # Add to global list of "all gaps".
#         global_gaps_all.extend(gaps_all_chrom)

#         # ----- Compute "F/R gaps" on this chromosome -----
#         # Sort sites by position, keeping strand info.
#         sites_sorted = sorted(all_sites_chrom, key=lambda x: x["pos"])

#         # Compute F/R or R/F gaps on this chromosome.
#         gaps_fr_chrom = compute_fr_gaps(
#             sites_sorted,
#             L,
#             circular=circular,
#         )

#         # Add to global list of F/R gaps.
#         global_gaps_fr.extend(gaps_fr_chrom)

#     # Build stats dictionary with global metrics.
#     stats = {
#         "n_primers": len(primers),
#         "n_chromosomes": n_chromosomes,
#         "n_chrom_with_sites": n_chrom_with_sites,
#         "n_sites": total_sites,
#         "n_gaps_all": len(global_gaps_all),
#         "mean_gap_all": safe_mean(global_gaps_all),
#         "max_gap_all": max(global_gaps_all) if global_gaps_all else math.nan,
#         "n_gaps_fr": len(global_gaps_fr),
#         "mean_gap_fr": safe_mean(global_gaps_fr),
#         "max_gap_fr": max(global_gaps_fr) if global_gaps_fr else math.nan,
#         "gaps_all": global_gaps_all,
#         "gaps_fr": global_gaps_fr,
#     }
#     return stats


def analyze_primer_set(
    primers,
    genome_dict,
    circular: bool,
    within_threshold: int,
    set_index: int | None = None,
    outdir: Path | None = None,
):
    """
    Analyze one primer set over a (possibly multi-chromosome) genome.

    - genome_dict: {chrom_name: chrom_seq}
    - circular: treat each chromosome as circular if True, else linear.
    - set_index, outdir: optional; if given, dump all binding sites to
      set_{set_index:03d}_sites.csv in outdir.

    For each chromosome:
        - Find binding sites (forward + reverse) for each primer.
        - Compute gaps (all) and F→R-only gaps.

    We do *not* connect chromosomes with gaps; each chromosome is independent.

    Returns:
        stats dict with keys:
            n_primers, n_chromosomes, n_chrom_with_sites,
            n_sites, n_gaps_all, mean_gap_all, max_gap_all,
            n_gaps_fr, mean_gap_fr, max_gap_fr,
            gaps_all, gaps_fr
    """
    # Total number of binding sites across all chromosomes.
    total_sites = 0

    # Total number of chromosomes in the genome.
    n_chromosomes = len(genome_dict)

    # Number of chromosomes that have at least one binding site.
    n_chrom_with_sites = 0

    # Global lists of gaps pooled across all chromosomes.
    global_gaps_all = []
    global_gaps_fr = []

    # SWGA2-style forward→reverse gaps (any F→R pair within 'within_threshold').
    # We also track per-chromosome max F→R gaps so we can average them like swga2.
    global_gaps_fr_within = []
    per_chrom_max_fr_within = []
    # global_gaps_fr_within: all F→R gaps within threshold across the genome.
    # per_chrom_max_fr_within: one max per chromosome.

    # Global list of all binding sites (for optional CSV dump).
    # Each entry: {"chrom": name, "pos": pos, "strand": "+/-", "primer": seq}
    all_sites_global = []

    # Iterate over chromosomes in the genome.
    for chrom, seq in genome_dict.items():
        L = len(seq)

        # Collect all binding sites for this primer set in this chromosome.
        all_sites_chrom = []
        for pr in primers:
            # find_sites returns [{"pos": int, "strand": "+/-"}, ...]
            for s in find_sites(seq, pr):
                # Enrich each site with the primer sequence and chromosome name.
                all_sites_chrom.append(
                    {
                        "chrom": chrom,
                        "pos": s["pos"],
                        "strand": s["strand"],
                        "primer": pr,
                    }
                )

        # If no sites on this chromosome, skip it.
        if not all_sites_chrom:
            continue

        # Record that this chromosome has sites.
        n_chrom_with_sites += 1
        total_sites += len(all_sites_chrom)

        # Add to global site list (for optional CSV later).
        all_sites_global.extend(all_sites_chrom)

        # ----- Compute "all gaps" on this chromosome -----
        positions = sorted(s["pos"] for s in all_sites_chrom)
        gaps_all_chrom = compute_gaps_from_positions(positions, L, circular=circular)
        global_gaps_all.extend(gaps_all_chrom)

        # ----- Compute F→R gaps on this chromosome -----
        sites_sorted = sorted(all_sites_chrom, key=lambda x: x["pos"])
        gaps_fr_chrom = compute_fr_gaps(sites_sorted, L, circular=circular)
        global_gaps_fr.extend(gaps_fr_chrom)

        # ----- SWGA2-style F→R gaps on this chromosome (any F→R within threshold) -----
        positions_fwd = [s["pos"] for s in all_sites_chrom if s["strand"] == "+"]
        positions_rev = [s["pos"] for s in all_sites_chrom if s["strand"] == "-"]

        gaps_fr_within_chrom = compute_fr_gaps_alternating(
            positions_fwd,
            positions_rev,
            threshold=within_threshold,
        )

        # Pool all F→R-within gaps across chromosomes.
        global_gaps_fr_within.extend(gaps_fr_within_chrom)

        # Track the per-chromosome max gap (to later average like swga2).
        if gaps_fr_within_chrom:
            per_chrom_max_fr_within.append(max(gaps_fr_within_chrom))

        # with this we have
        # - global_gaps_fr_within: the global list of all F→R distances (within threshold).
        # - per_chrom_max_fr_within: one max gap per chromosome.
        # ------------------------------------------------------------

    # ----- Optional: dump all binding sites for this set to CSV -----
    if set_index is not None and outdir is not None and all_sites_global:
        import csv

        sites_csv = outdir / f"set_{set_index:03d}_sites.csv"
        # Sort for nicer viewing: by chromosome, then position.
        all_sites_global_sorted = sorted(
            all_sites_global, key=lambda s: (s["chrom"], s["pos"])
        )
        with sites_csv.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["chrom", "pos", "strand", "primer"])
            for s in all_sites_global_sorted:
                writer.writerow([s["chrom"], s["pos"], s["strand"], s["primer"]])
    # ----------------------------------------------------------------

    # Build stats dictionary with global metrics.

    # SWGA2-style F→R summary:
    # - global_gaps_fr_within: all F→R gaps within threshold across chromosomes.
    # - per_chrom_max_fr_within: per-chromosome max F→R gaps (to average like swga2).
    if per_chrom_max_fr_within:
        fg_max_forward_reverse_mean_max_per_chr = sum(per_chrom_max_fr_within) / len(
            per_chrom_max_fr_within
        )
    else:
        fg_max_forward_reverse_mean_max_per_chr = math.nan

    if global_gaps_fr_within:
        fg_max_forward_reverse_global = max(global_gaps_fr_within)
    else:
        fg_max_forward_reverse_global = math.nan

    stats = {
        "n_primers": len(primers),
        "n_chromosomes": n_chromosomes,
        "n_chrom_with_sites": n_chrom_with_sites,
        "n_sites": total_sites,
        "n_gaps_all": len(global_gaps_all),
        "mean_gap_all": safe_mean(global_gaps_all),
        "max_gap_all": max(global_gaps_all) if global_gaps_all else math.nan,
        "n_gaps_fr": len(global_gaps_fr),
        "mean_gap_fr": safe_mean(global_gaps_fr),
        "max_gap_fr": max(global_gaps_fr) if global_gaps_fr else math.nan,
        # SWGA2-style F→R gaps:
        "n_gaps_fr_within": len(global_gaps_fr_within),
        "mean_gap_fr_within": safe_mean(global_gaps_fr_within),
        "fg_max_forward_reverse_mean_max_per_chr": fg_max_forward_reverse_mean_max_per_chr,
        "fg_max_forward_reverse_global": fg_max_forward_reverse_global,
        "gaps_all": global_gaps_all,
        "gaps_fr": global_gaps_fr,
    }
    return stats


# ===========================
# Plotting
# ===========================


def plot_for_set(idx, score, stats, outdir: Path):
    """
    Generate plots for one primer set:

    - Histogram of all gaps (global across chromosomes)
    - Histogram comparing all gaps vs F/R-only gaps
    - Sorted gaps plot (all gaps, global)

    Files are saved in 'outdir' as PNGs with set index in the filename.
    """
    # Extract global gaps lists from stats.
    gaps_all = stats["gaps_all"]
    gaps_fr = stats["gaps_fr"]

    # If there are no gaps (no sites at all), we do not plot anything.
    if not gaps_all:
        print(f"[Set {idx}] No gaps to plot (no binding sites in any chromosome).")
        return

    # # ----- 1) Histogram of all gaps -----
    # plt.figure()
    # plt.hist(gaps_all, bins=20)
    # plt.xlabel("Gap length (bp)")
    # plt.ylabel("Count")
    # plt.title(f"Set {idx}: all gaps (global; score={score})")
    # plt.tight_layout()
    # # Example filename: set_001_hist_all.png
    # plt.savefig(outdir / f"set_{idx:03d}_hist_all.png")
    # plt.close()

    # ----- 2) Histogram: all gaps vs F/R-only gaps -----
    plt.figure()
    # Histogram of all gaps
    plt.hist(gaps_all, bins=20, alpha=0.5, label="All gaps")
    # If we have F/R gaps, overlay them.
    if gaps_fr:
        plt.hist(gaps_fr, bins=20, alpha=0.5, label="F/R gaps")
    plt.xlabel("Gap length (bp)")
    plt.ylabel("Count")
    plt.legend()
    plt.title(f"Set {idx}: all vs F/R gaps (global; score={score})")
    plt.tight_layout()
    plt.savefig(outdir / f"set_{idx:03d}_hist_all_vs_fr.png")
    plt.close()

    # # ----- 3) Sorted gaps (all gaps) -----
    # # Sort gaps from smallest to largest.
    # gaps_sorted = sorted(gaps_all)

    # plt.figure()
    # # x-axis is the rank (1, 2, ..., n), y-axis is the sorted gap length.
    # plt.plot(range(1, len(gaps_sorted) + 1), gaps_sorted, marker="o")
    # plt.xlabel("Gap rank (small → large)")
    # plt.ylabel("Gap length (bp)")
    # plt.title(f"Set {idx:03d}: sorted gaps (global; score={score})")
    # plt.tight_layout()
    # plt.savefig(outdir / f"set_{idx:03d}_sorted_gaps.png")
    # plt.close()


# ===========================
# Main entry point
# ===========================


def main():
    """
    Main function:

    - Parse arguments.
    - Load genome (multi-FASTA) into a dict.
    - Load primer sets from input file (step4 or generic CSV).
    - For each primer set:
        * Analyze gaps.
        * Collect summary stats.
        * Generate global plots.
    - Save summary statistics to a CSV file.
    """
    # Parse command-line options.
    args = parse_args()

    # Ensure the output directory exists (create it if necessary).
    outdir: Path = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    # Load the genome as a dict {seq_name: seq}.
    genome_dict = load_fasta_dict(args.input_genome)

    # Print some basic info about the genome.
    print(
        f"Loaded genome from {args.input_genome} "
        f"({len(genome_dict)} sequences: "
        f"{', '.join(list(genome_dict.keys())[:5])}"
        f"{'...' if len(genome_dict) > 5 else ''})"
    )
    print(f"Circular sequences: {args.circular}")

    # Decide how to parse the primer sets based on --step4-format flag.
    if args.step4_format:
        primer_sets = load_primer_sets_step4(args.input_file)
    else:
        primer_sets = load_primer_sets_generic_csv(args.input_file)

    # Report how many primer sets were loaded.
    print(f"Loaded {len(primer_sets)} primer sets from {args.input_file}")

    results = []  # Will store one summary row per primer set.

    # Iterate over primer sets with an index starting at 1.
    for idx, (primers, score) in enumerate(primer_sets, start=1):
        print(f"\nAnalyzing set {idx} (score={score}, {len(primers)} primers)...")

        # Compute binding sites and gap statistics for this set over the genome.
        stats = analyze_primer_set(
            primers=primers,
            genome_dict=genome_dict,
            circular=args.circular,
            within_threshold=args.within_threshold,
            set_index=idx,
            outdir=outdir,
        )

        # ---- DEBUG: inspect F→R gaps for a specific set ----
        if idx == 1:  # change 1 to 2,3,... to debug another set
            gaps_fr = stats["gaps_fr"]
            print(f"[DEBUG] Set {idx}: F→R gaps ({len(gaps_fr)} total):")
            print(gaps_fr)
            if gaps_fr:
                print(f"[DEBUG] max(gaps_fr) = {max(gaps_fr)}")
                print(f"[DEBUG] stats['max_gap_fr'] = {stats['max_gap_fr']}")
        # ---- END DEBUG ----

        # Build a "flat" summary row (no lists, only scalar numbers).
        stats_row = {
            "set_index": idx,
            "score": score,
            "n_primers": stats["n_primers"],
            "n_chromosomes": stats["n_chromosomes"],
            "n_chrom_with_sites": stats["n_chrom_with_sites"],
            "n_sites": stats["n_sites"],
            "n_gaps_all": stats["n_gaps_all"],
            "mean_gap_all": stats["mean_gap_all"],
            "max_gap_all": stats["max_gap_all"],
            "n_gaps_fr": stats["n_gaps_fr"],
            "mean_gap_fr": stats["mean_gap_fr"],
            "max_gap_fr": stats["max_gap_fr"],
            "n_gaps_fr_within": stats["n_gaps_fr_within"],
            "mean_gap_fr_within": stats["mean_gap_fr_within"],
            "fg_max_forward_reverse_mean_max_per_chr": stats[
                "fg_max_forward_reverse_mean_max_per_chr"
            ],
            "fg_max_forward_reverse_global": stats["fg_max_forward_reverse_global"],
        }

        results.append(stats_row)

        # Generate plots for this primer set (global across chromosomes).
        plot_for_set(idx, score, stats, outdir)

    # Convert the list of summary rows to a DataFrame.
    df_res = pd.DataFrame(results)

    # Output CSV path (inside the output directory).
    out_csv = outdir / "gap_stats_per_set.csv"

    # Save summary table to CSV (no index column).
    df_res.to_csv(out_csv, index=False)

    # Final messages to the user.
    print(f"\nSummary saved to {out_csv}")
    print(f"Plots saved to {outdir}/set_XXX_*.png")


# Run main() when the script is executed as a program.
if __name__ == "__main__":
    main()
