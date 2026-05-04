import re
import pandas as pd

# =========================
# 1. STRUCTURE PARSING
# =========================

def extract_ss_cons(stockholm_file):
    """
    Extract SS_cons line from Stockholm file.
    """
    with open(stockholm_file) as f:
        for line in f:
            if "SS_cons" in line:
                return line.split()[-1]
    return None


def dotbracket_to_pairs(structure):
    """
    Convert dot-bracket to base pairs.
    """
    stack = []
    pairs = set()

    for i, c in enumerate(structure):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            pairs.add((j, i))

    return pairs


# =========================
# 2. COV FILE PARSING (ROBUST)
# =========================

def parse_cov_generic(file_path):
    """
    Flexible parser for .cov-like files.
    Does NOT assume strict format.
    Extracts ANY lines containing two indices + numeric score/E-value.
    """

    rows = []

    with open(file_path) as f:
        for line in f:
            line = line.strip()

            # skip headers/comments
            if not line or line.startswith("#"):
                continue

            parts = line.split()

            # try to detect numeric structure
            nums = []
            for p in parts:
                try:
                    nums.append(float(p))
                except:
                    pass

            # we expect at least i, j
            if len(nums) >= 2:
                i, j = int(nums[0]), int(nums[1])

                score = nums[2] if len(nums) > 2 else None
                evalue = nums[3] if len(nums) > 3 else None

                rows.append((i, j, score, evalue))

    return pd.DataFrame(rows, columns=["i", "j", "score", "evalue"])


# =========================
# 3. STRUCTURE COMPARISON
# =========================

def compare_structures(pred_pairs, ref_pairs):
    pred = set(pred_pairs)
    ref = set(ref_pairs)

    tp = len(pred & ref)
    fp = len(pred - ref)
    fn = len(ref - pred)

    precision = tp / (tp + fp) if (tp + fp) else 0
    recall = tp / (tp + fn) if (tp + fn) else 0

    return {
        "TP": tp,
        "FP": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall
    }


# =========================
# 4. COVARIATION SUPPORT
# =========================

def covariation_support(struct_pairs, cov_df, evalue_threshold=0.05):
    """
    Fraction of predicted pairs supported by covariation.
    """

    if cov_df.empty:
        return 0.0

    cov_pairs = set(
        (row.i, row.j)
        for row in cov_df.itertuples()
        if row.evalue is None or row.evalue <= evalue_threshold
    )

    supported = len(struct_pairs & cov_pairs)

    return supported / len(struct_pairs) if struct_pairs else 0


# =========================
# 5. FULL METHOD EVALUATION
# =========================

def evaluate_method(struct_file, cov_file, ref_structure_file):
    """
    Main evaluation function for CaCoFold or RNAalifold output.
    """

    # --- predicted structure ---
    pred_structure = extract_ss_cons(struct_file)
    if pred_structure is None:
        raise ValueError(f"No SS_cons found in {struct_file}")

    pred_pairs = dotbracket_to_pairs(pred_structure)

    # --- covariation ---
    cov_df = parse_cov_generic(cov_file)

    # --- reference structure ---
    ref_structure = extract_ss_cons(ref_structure_file)
    if ref_structure is None:
        raise ValueError(f"No SS_cons in reference {ref_structure_file}")

    ref_pairs = dotbracket_to_pairs(ref_structure)

    # --- metrics ---
    struct_metrics = compare_structures(pred_pairs, ref_pairs)
    cov_support = covariation_support(pred_pairs, cov_df)

    return {
        **struct_metrics,
        "covariation_support": cov_support,
        "n_pairs": len(pred_pairs)
    }


# =========================
# 6. RUN EXAMPLE
# =========================

# FILE PATHS (CHANGE THESE)
struct_file = "RF00005_tRNA.cacofold.sto"
cov_file = "RF00005_tRNA.cacofold.cov"
ref_file = "RF00005_tRNA.R2R.sto"

results = evaluate_method(struct_file, cov_file, ref_file)

print(results)