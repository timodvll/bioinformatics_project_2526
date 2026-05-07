from Bio import AlignIO
import random


#RQ 1:

def extract_consensus_structure(sto_file):
    with open(sto_file) as f:
        for line in f:
            if line.startswith("#=GC SS_cons"):
                return line.strip().split()[2]
    return None

def extract_cacofold_structure(sto_file):
    with open(sto_file) as f:
        for line in f:
            if "SS_cons" in line:
                parts = line.strip().split()
                if len(parts) >= 3:
                    return parts[2]
    return None


def dotbracket_to_pairs(structure):
    stack = {}
    pairs = set()

    opening = "([{<"
    closing = ")]}>"

    for i, c in enumerate(structure, start=1):
        if c in opening:
            if c not in stack:
                stack[c] = []
            stack[c].append(i)

        elif c in closing:
            open_c = opening[closing.index(c)]
            if open_c in stack and stack[open_c]:
                j = stack[open_c].pop()
                pairs.add((j, i))

    return pairs

def count_pseudoknots(pairs):
    pairs = sorted(pairs)
    pk = 0

    for i1, j1 in pairs:
        for i2, j2 in pairs:
            if i1 < i2 < j1 < j2:
                pk += 1

    return pk // 2  # avoid double counting


def compare_structures(consensus, cacofold):
    c_pairs = dotbracket_to_pairs(consensus)
    k_pairs = dotbracket_to_pairs(cacofold)

    c_knots = count_pseudoknots(c_pairs)
    k_knots = count_pseudoknots(k_pairs)

    c_helices = extract_helices(c_pairs)
    c_helices_stats = helix_stats(c_helices)

    k_helices = extract_helices(k_pairs)
    k_helices_stats = helix_stats(k_helices)

    # ---- core metrics ----
    tp_pairs = c_pairs & k_pairs
    fp_pairs = k_pairs - c_pairs
    fn_pairs = c_pairs - k_pairs

    tp = len(tp_pairs)
    fp = len(fp_pairs)
    fn = len(fn_pairs)

    precision = tp / (tp + fp) if (tp + fp) else 0
    recall = tp / (tp + fn) if (tp + fn) else 0

    return {
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall,
        "c_knots": c_knots,
        "k_knots": k_knots,
        "c_helices": c_helices_stats,
        "k_helices": k_helices_stats,

        "fp_pairs": sorted(list(fp_pairs)),
        "tp_pairs": sorted(list(tp_pairs)),
        "fn_pairs": sorted(list(fn_pairs))
    }

def helix_stats(helices):
    lengths = [len(h) for h in helices]

    return {
        "num_helices": len(helices),
        "avg_length": sum(lengths) / len(lengths) if lengths else 0,
        "max_length": max(lengths) if lengths else 0
    }

def extract_helices(pairs):
    pairs = sorted(pairs)
    used = set()
    helices = []

    for i, (a, b) in enumerate(pairs):
        if (a, b) in used:
            continue

        helix = [(a, b)]
        used.add((a, b))

        # extend helix
        next_a, next_b = a + 1, b - 1
        while (next_a, next_b) in pairs:
            helix.append((next_a, next_b))
            used.add((next_a, next_b))
            next_a += 1
            next_b -= 1

        helices.append(helix)

    return helices

def get_fp_pairs(result_dict):
    """
    Extract false positive base pairs from compare_structures output.
    """
    return set(result_dict["fp_pairs"])

def map_fp_to_cov(fp_pairs, cov_df, window=3):
    """
    Match FP pairs to cov pairs allowing small coordinate shifts.
    """

    cov_list = list(zip(cov_df["left_pos"], cov_df["right_pos"]))

    results = []

    for i, j in fp_pairs:

        found = None

        for ci, cj in cov_list:

            if abs(i - ci) <= window and abs(j - cj) <= window:
                found = cov_df[
                    (cov_df["left_pos"] == ci) &
                    (cov_df["right_pos"] == cj)
                ].iloc[0]
                break

        if found is not None:
            results.append({
                "pair": (i, j),
                "evalue": found["evalue"],
                "supported": found["evalue"] < 0.05
            })
        else:
            results.append({
                "pair": (i, j),
                "evalue": None,
                "supported": False
            })

    return results

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


# --------------------------------------------------
# Parsing
# --------------------------------------------------
def parse_cacofold_cov(filepath):
    rows = []

    with open(filepath, "r") as f:
        for line in f:

            if (
                not line.strip()
                or line.startswith("#")
                or line.startswith("-")
            ):
                continue

            parts = line.rstrip("\n").split()

            if len(parts) == 8:
                in_cacofold = parts[0]
                in_given = ""
                offset = 1

            elif len(parts) >= 9:
                in_cacofold = parts[0]
                in_given = parts[1]
                offset = 2

            else:
                continue

            rows.append({
                "in_cacofold": in_cacofold,
                "in_given": in_given,
                "left_pos": int(parts[offset]),
                "right_pos": int(parts[offset + 1]),
                "score": float(parts[offset + 2]),
                "evalue": float(parts[offset + 3]),
                "pvalue": float(parts[offset + 4]),
                "substitutions": int(parts[offset + 5]),
                "power": float(parts[offset + 6]),
            })

    df = pd.DataFrame(rows)

    df["left_pos"] = df["left_pos"].astype(int)
    df["right_pos"] = df["right_pos"].astype(int)

    return df

# --------------------------------------------------
# Extra pair extraction
# --------------------------------------------------
def get_extra_pairs(df):
    """
    Strict extra pairs:
    predicted by CaCoFold but absent from consensus.
    """

    return df[
        (df["in_cacofold"] == "*") &
        (df["in_given"] == "")
    ].copy()


# --------------------------------------------------
# Biological support classification
# --------------------------------------------------

def classify_support(df):
    """
    Add support labels based on E-value
    """

    df = df.copy()

    def label(ev):
        if ev < 0.01:
            return "highly_significant"
        elif ev < 0.05:
            return "significant"
        else:
            return "unsupported"

    df["support"] = df["evalue"].apply(label)

    return df


# --------------------------------------------------
# Summary metrics
# --------------------------------------------------

def summarize_extra_pairs(df):
    """
    Compute statistics for extra pairs
    """

    total = len(df)

    high = (df["support"] == "highly_significant").sum()
    sig = (df["support"] == "significant").sum()
    unsup = (df["support"] == "unsupported").sum()

    return {
        "total_extra_pairs": total,
        "highly_significant": high,
        "significant": sig,
        "unsupported": unsup,
        "supported_fraction": (high + sig) / total if total else 0
    }


# --------------------------------------------------
# Plotting
# --------------------------------------------------

def plot_support_distribution(df):
    counts = df["support"].value_counts()

    counts.plot(kind="bar")
    plt.ylabel("Count")
    plt.title("Covariation Support of Extra CaCoFold Pairs")
    plt.show()


def plot_evalue_histogram(df):
    plt.hist(df["evalue"], bins=40)
    plt.xscale("log")
    plt.xlabel("E-value (log scale)")
    plt.ylabel("Count")
    plt.title("Distribution of Extra Pair E-values")
    plt.show()


# --------------------------------------------------
# Full pipeline
# --------------------------------------------------

def analyze_cacofold_cov(filepath):
    df = parse_cacofold_cov(filepath)
    extra = get_extra_pairs(df)
    extra = classify_support(extra)
    summary = summarize_extra_pairs(extra)

    return df, extra, summary


#RQ 2:
from pathlib import Path


from pathlib import Path
import random

def subsample_interleaved_stockholm(input_file, fraction, seed):
    random.seed(seed)
    input_file = Path(input_file)

    with open(input_file) as f:
        lines = f.readlines()

    # collect sequence IDs
    seq_ids = []

    for line in lines:
        if line.startswith("#") or line.strip() == "" or line.strip() == "//":
            continue

        seq_id = line.split()[0]

        if seq_id not in seq_ids:
            seq_ids.append(seq_id)

    n_keep = max(1, int(len(seq_ids) * fraction))
    keep_ids = set(random.sample(seq_ids, n_keep))

    print(f"Keeping {n_keep}/{len(seq_ids)} sequences")

    # -----------------------------
    # ONLY CHANGE: output structure
    # -----------------------------
    frac_dir = f"{int(fraction * 100)}"
    out_dir = Path("RQ2/05") / frac_dir / str(seed)
    out_dir.mkdir(parents=True, exist_ok=True)

    outname = input_file.stem + f"_{int(fraction * 100)}_{seed}.sto"
    outpath = out_dir / outname

    with open(outpath, "w") as out:
        for line in lines:

            if (
                line.startswith("#=GC")
                or line.startswith("# STOCKHOLM")
                or line.strip() == "//"
                or line.strip() == ""
            ):
                out.write(line)

            elif line.startswith("#=GR"):
                seq_id = line.split()[1]
                if seq_id in keep_ids:
                    out.write(line)

            elif line.startswith("#"):
                out.write(line)

            else:
                seq_id = line.split()[0]
                if seq_id in keep_ids:
                    out.write(line)

    print(f"Wrote {outpath}")
    return str(outpath)


# RQ3
from pathlib import Path


def local_swap(seq):
    """
    Perform one local ±1 alignment swap.
    """
    seq = list(seq)

    valid_positions = [
        i for i in range(1, len(seq))
        if seq[i] != "-" and seq[i-1] != "-"
    ]

    if valid_positions:
        i = random.choice(valid_positions)
        seq[i], seq[i-1] = seq[i-1], seq[i]

    return "".join(seq)


def perturb_sequence(seq, n_shifts=5):
    """
    Apply multiple local swaps to a sequence.
    """
    for _ in range(n_shifts):
        seq = local_swap(seq)
    return seq


def perturb_stockholm_local(input_file, noise_percent, seed):
    random.seed(seed)

    input_file = Path(input_file)

    with open(input_file) as f:
        lines = f.readlines()

    # Collect unique sequence IDs
    seq_ids = []

    for line in lines:
        if line.startswith("#") or line.strip() == "" or line.strip() == "//":
            continue

        seq_id = line.split()[0]

        if seq_id not in seq_ids:
            seq_ids.append(seq_id)

    n_perturb = max(1, int(len(seq_ids) * noise_percent / 100))
    perturbed_ids = set(random.sample(seq_ids, n_perturb))

    outpath = Path(f"RQ3_2/05/{noise_percent}/{seed}/RF00005_{noise_percent}_{seed}.sto")

    with open(outpath, "w") as out:
        for line in lines:

            if line.startswith("#") or line.strip() == "" or line.strip() == "//":
                out.write(line)

            else:
                parts = line.rstrip("\n").split()

                seq_id = parts[0]
                seq = parts[-1]

                if seq_id in perturbed_ids:
                    seq = perturb_sequence(seq, n_shifts=50)

                prefix = line[:line.find(parts[-1])]
                out.write(prefix + seq + "\n")

    print(
        f"Created {outpath} "
        f"({n_perturb}/{len(seq_ids)} sequences perturbed, 5 shifts each)"
    )


# RQ4
