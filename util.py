from Bio import AlignIO


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

    for i, c in enumerate(structure):
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


    tp = len(c_pairs & k_pairs)
    fp = len(k_pairs - c_pairs)
    fn = len(c_pairs - k_pairs)

    precision = tp / (tp + fp) if (tp + fp) else 0
    recall = tp / (tp + fn) if (tp + fn) else 0

    return {
        "TP": tp,
        "FP (extra CaCoFold)": fp,
        "FN": fn,
        "precision": precision,
        "recall": recall,
        "c_knots": c_knots,
        "k_knots": k_knots,
        "c_helices": c_helices_stats,
        "k_helices": k_helices_stats
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


#RQ 2:
import random
from pathlib import Path


def subsample_interleaved_stockholm(input_file, fraction=0.5, seed=42):
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

    outname = input_file.stem + f"_{int(fraction * 100)}.sto"

    with open(outname, "w") as out:
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

    print(f"Wrote {outname}")
    return outname

# RQ3
import random
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

    outpath = Path(f"RQ3/{noise_percent}/{seed}/RF00162_{noise_percent}_{seed}.sto")

    with open(outpath, "w") as out:
        for line in lines:

            if line.startswith("#") or line.strip() == "" or line.strip() == "//":
                out.write(line)

            else:
                parts = line.rstrip("\n").split()

                seq_id = parts[0]
                seq = parts[-1]

                if seq_id in perturbed_ids:
                    seq = perturb_sequence(seq, n_shifts=5)

                prefix = line[:line.find(parts[-1])]
                out.write(prefix + seq + "\n")

    print(
        f"Created {outpath} "
        f"({n_perturb}/{len(seq_ids)} sequences perturbed, 5 shifts each)"
    )

    import random
from pathlib import Path


def shift_block(seq_list, start, end, direction):
    """
    Shift a contiguous block by ±1 position.
    """
    if direction == "left" and start > 0:
        block = seq_list[start:end]
        seq_list[start-1:start-1] = block
        del seq_list[end]

    elif direction == "right" and end < len(seq_list):
        block = seq_list[start:end]
        del seq_list[start:end]
        seq_list[start+1:start+1] = block

    return seq_list


def apply_block_misalignment(seq, n_blocks=2):
    """
    Apply block-level alignment perturbations.
    Each block: length 3–8, shift ±1.
    """
    seq = list(seq)

    for _ in range(n_blocks):

        L = random.randint(3, 8)

        if len(seq) <= L + 2:
            continue

        start = random.randint(1, len(seq) - L - 1)
        end = start + L

        direction = random.choice(["left", "right"])

        seq = shift_block(seq, start, end, direction)

    return "".join(seq)


def perturb_stockholm_block(input_file, noise_percent, seed):
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

    n_perturb = max(1, int(len(seq_ids) * noise_percent / 100))
    perturbed_ids = set(random.sample(seq_ids, n_perturb))

    outpath = Path(
        f"RQ3_2/{noise_percent}/{seed}/RF00162_{noise_percent}_{seed}.sto"
    )

    with open(outpath, "w") as out:
        for line in lines:

            if line.startswith("#") or line.strip() == "" or line.strip() == "//":
                out.write(line)

            else:
                parts = line.rstrip("\n").split()
                seq_id = parts[0]
                seq = parts[-1]

                if seq_id in perturbed_ids:
                    seq = apply_block_misalignment_safe(seq, n_blocks=2)

                prefix = line[:line.find(parts[-1])]
                out.write(prefix + seq + "\n")

    print(
        f"Created {outpath} "
        f"({n_perturb}/{len(seq_ids)} sequences, block misalignment)"
    )

import random


def apply_block_misalignment_safe(seq, n_blocks=2):
    """
    Safe block misalignment that preserves alignment length.
    """
    seq = list(seq)

    for _ in range(n_blocks):

        L = random.randint(3, 8)

        if len(seq) <= L + 2:
            continue

        start = random.randint(1, len(seq) - L - 1)
        end = start + L

        direction = random.choice(["left", "right"])

        block = seq[start:end]

        # replace original block with gaps
        seq[start:end] = ["-"] * L

        if direction == "left" and start > 0:
            seq[start-1:start-1] = block[:-1]

        elif direction == "right" and end < len(seq):
            seq[end:end] = block[1:]

    return "".join(seq)