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