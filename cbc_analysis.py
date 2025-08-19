import re
import argparse
import sys
import csv

# =========================
# Normalization helpers
# =========================


def normalize_base(b):
    """Normalize base for ITS2: T->U, keep gaps and other chars as-is."""
    if not b:
        return b
    b = b.upper()
    return 'U' if b == 'T' else b


def normalize_arrow(s):
    """Normalize change arrows to a single '->' representation."""
    if s is None:
        return ""
    return re.sub(r"(==>|-->|→|⇒|—>|->)", "->", s)


# =========================
# Parsers
# =========================


def parse_xfasta(path):
    """
    Parse XFasta: >name, aligned SEQ, aligned STRUCT (dot-bracket)
    """
    taxa = {}
    with open(path, "r") as f:
        lines = [l.rstrip("\n") for l in f]
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            name = lines[i][1:].strip().replace(" ", "_")
            if i + 2 >= len(lines):
                raise ValueError(f"{name} has incomplete seq/structure lines.")
            seq = lines[i+1].strip()
            struct = lines[i+2].strip()
            if len(seq) != len(struct):
                raise ValueError(f"[{name}] sequence vs structure length mismatch: {len(seq)} vs {len(struct)}")
            taxa[name] = (seq, struct)
            i += 3
        else:
            i += 1
    return taxa


def parse_changes_with_events(changes_path):
    """
    Parse PAUP 'Apomorphy lists' block (from 'describetrees / apolist=yes')
    into: pos (1-based alignment column) -> [event lines...]
    Each event line will include the branch header for context.
    """
    with open(changes_path, "r") as f:
        lines = [ln.rstrip("\n") for ln in f]

    # Find start of "Apomorphy lists:"
    start_idx = None
    for i, ln in enumerate(lines):
        if re.search(r"^\s*Apomorphy lists:\s*$", ln):
            start_idx = i + 1
            break
    if start_idx is None:
        raise ValueError("Could not find 'Apomorphy lists:' section in PAUP output. Make sure you used 'describetrees / apolist=yes'.")

    # Patterns
    # Example header: "node_6 --> node_5" or "node_5 --> tetrix bolivari"
    header_pat = re.compile(r"^\s*(node_\d+\s+-->\s+(?:node_\d+|[A-Za-z_ ]+\(?\d*\)?))\s*$")
    # Example entry lines (columns): "   10     1   1.000  C ==> A" or "... A --> C"
    entry_pat = re.compile(r"^\s*(\d+)\s+\d+\s+[\d.]+\s+([ACGTU\-])\s*(?:==>|-->|→|⇒|->)\s*([ACGTU\-])\s*$")

    header_entry_pat = re.compile(
        r"^\s*(node_\d+\s+-->\s+(?:node_\d+|[A-Za-z_ ]+\(?\d*\)?))\s+"
        r"(\d+)\s+\d+\s+[\d.]+\s+([ACGTU\-])\s*(?:==>|-->|→|⇒|->)\s*([ACGTU\-])\s*$"
    )

    pos2events = {}
    current_header = None

    # Skip down to the table body (under the dashed line)
    i = start_idx
    while i < len(lines) and not re.search(r"^-{3,}", lines[i]):
        mhead = header_pat.match(lines[i])
        if mhead:
            current_header = mhead.group(1)
        i += 1

    # Walk through the sections
    for j in range(i, len(lines)):
        ln = lines[j]
        if not ln.strip():
            continue
        mheadrow = header_entry_pat.match(ln)
        if mheadrow:
            current_header = mheadrow.group(1)
            pos = int(mheadrow.group(2))
            b_from = mheadrow.group(3)
            b_to = mheadrow.group(4)
            ev = f"{current_header}: {b_from} -> {b_to}"
            pos2events.setdefault(pos, []).append(ev)
            continue
        mhead = header_pat.match(ln)
        if mhead:
            current_header = mhead.group(1)
            continue
        mrow = entry_pat.match(ln)
        if mrow and current_header:
            pos = int(mrow.group(1))
            b_from = mrow.group(2)
            b_to   = mrow.group(3)
            ev = f"{current_header}: {b_from} -> {b_to}"
            pos2events.setdefault(pos, []).append(ev)
            continue
        if re.search(r"^-{3,}", ln):
            continue
    return pos2events

# =========================
# Structure helpers
# =========================


def mask_structure_by_gaps(seq, struct):
    """Replace '(' or ')' with '.' in columns where there is a gap in the sequence."""
    out = list(struct)
    for i, ch in enumerate(out):
        if seq[i] == '-' and ch in '()':
            out[i] = '.'
    return ''.join(out)


def find_pairs(struct_masked):
    """Return base-pair mapping from dot-bracket string (0-based indexes)."""
    stack, pairs = [], {}
    for i, ch in enumerate(struct_masked):
        if ch == '(':
            stack.append(i)
        elif ch == ')':
            if stack:
                j = stack.pop()
                pairs[i] = j
                pairs[j] = i
    return pairs


def build_ungapped_prefix(seq):
    """Return prefix sums: prefix[i] = number of real bases up to and including index i (0-based)."""
    prefix = []
    count = 0
    for ch in seq:
        if ch.upper() in ("A", "U", "G", "C", "T"):
            count += 1
        prefix.append(count)
    return prefix


def ungapped_at(prefix, idx0):
    """Get 1-based ungapped position using precomputed prefix array."""
    if idx0 < 0 or idx0 >= len(prefix):
        return ""
    return prefix[idx0]


def ungapped_pos_1based(seq, aln_idx0, prefix=None):
    """Convert alignment index (0-based) to ungapped 1-based position (O(1) if prefix given)."""
    if prefix is not None:
        return ungapped_at(prefix, aln_idx0)
    # Fallback (should not be used in hot paths)
    count = 0
    for i, ch in enumerate(seq):
        if ch.upper() in ("A", "U", "G", "C", "T"):
            count += 1
        if i == aln_idx0:
            return count
    return ""


def valid_pair(b1, b2):
    """Check if two bases form a valid canonical or GU pair."""
    s = {normalize_base(b1), normalize_base(b2)}
    return s in ({"A", "U"}, {"G", "C"}, {"G", "U"})


def classify_cbc_at_cached(aln_pos_1based, seq1, struct1, seq2, struct2,
                           mask1, mask2, pairs1, pairs2, ung1_prefix, ung2_prefix):
    """
    Cached version: uses precomputed mask/pairs and ungapped prefix arrays.
    """
    i = aln_pos_1based - 1
    base1, base2 = seq1[i], seq2[i]

    partner1 = pairs1.get(i)
    partner2 = pairs2.get(i)

    paired1 = (mask1[i] in '()') and (base1 != '-') and (partner1 is not None) and (seq1[partner1] != '-')
    paired2 = (mask2[i] in '()') and (base2 != '-') and (partner2 is not None) and (seq2[partner2] != '-')

    ung1 = ungapped_pos_1based(seq1, i, ung1_prefix)
    ung2 = ungapped_pos_1based(seq2, i, ung2_prefix)

    if not paired1 or not paired2:
        status = "Gap" if (base1 == '-' or base2 == '-') else "Unpaired"
        return {
            "s1_base": base1, "s1_struct": struct1[i], "s1_aln": aln_pos_1based,
            "s1_ung": ung1, "s1_pair_aln": "", "s1_pair_ung": "", "s1_pair_base": "",
            "s2_base": base2, "s2_struct": struct2[i], "s2_aln": aln_pos_1based,
            "s2_ung": ung2, "s2_pair_aln": "", "s2_pair_ung": "", "s2_pair_base": "",
            "type": status
        }

    partner1_base = seq1[partner1]
    partner2_base = seq2[partner2]

    pair1_bases = {normalize_base(base1), normalize_base(partner1_base)}
    pair2_bases = {normalize_base(base2), normalize_base(partner2_base)}

    pair1_aln = partner1 + 1
    pair2_aln = partner2 + 1
    pair1_ung = ungapped_pos_1based(seq1, partner1, ung1_prefix)
    pair2_ung = ungapped_pos_1based(seq2, partner2, ung2_prefix)

    if pair1_bases == pair2_bases:
        typ = "No Change"
    else:
        v1, v2 = valid_pair(*pair1_bases), valid_pair(*pair2_bases)
        if v1 and v2:
            changed_site = (normalize_base(base1) != normalize_base(base2))
            changed_partner = (normalize_base(partner1_base) != normalize_base(partner2_base))
            typ = "CBC" if (changed_site and changed_partner) else "hCBC"
        elif v1 or v2:
            typ = "hCBC"
        else:
            typ = "No Pair"

    return {
        "s1_base": base1, "s1_struct": struct1[i], "s1_aln": aln_pos_1based,
        "s1_ung": ung1, "s1_pair_aln": pair1_aln, "s1_pair_ung": pair1_ung, "s1_pair_base": partner1_base,
        "s2_base": base2, "s2_struct": struct2[i], "s2_aln": aln_pos_1based,
        "s2_ung": ung2, "s2_pair_aln": pair2_aln, "s2_pair_ung": pair2_ung, "s2_pair_base": partner2_base,
        "type": typ
    }


def classify_cbc_at(aln_pos_1based, seq1, struct1, seq2, struct2):
    """
    Classify CBC/hCBC/Gap/Unpaired at a given alignment position.
    Includes aligned and ungapped coordinates for both site and its partner.
    """
    i = aln_pos_1based - 1
    base1, base2 = seq1[i], seq2[i]

    mask1 = mask_structure_by_gaps(seq1, struct1)
    mask2 = mask_structure_by_gaps(seq2, struct2)
    pairs1 = find_pairs(mask1)
    pairs2 = find_pairs(mask2)
    partner1 = pairs1.get(i)
    partner2 = pairs2.get(i)

    paired1 = (mask1[i] in '()') and (base1 != '-') and (partner1 is not None) and (seq1[partner1] != '-')
    paired2 = (mask2[i] in '()') and (base2 != '-') and (partner2 is not None) and (seq2[partner2] != '-')

    ung1 = ungapped_pos_1based(seq1, i)
    ung2 = ungapped_pos_1based(seq2, i)

    if not paired1 or not paired2:
        status = "Gap" if (base1 == '-' or base2 == '-') else "Unpaired"
        return {
            "s1_base": base1, "s1_struct": struct1[i], "s1_aln": aln_pos_1based,
            "s1_ung": ung1, "s1_pair_aln": "", "s1_pair_ung": "", "s1_pair_base": "",
            "s2_base": base2, "s2_struct": struct2[i], "s2_aln": aln_pos_1based,
            "s2_ung": ung2, "s2_pair_aln": "", "s2_pair_ung": "", "s2_pair_base": "",
            "type": status
        }

    pair1_bases = {base1.upper(), seq1[partner1].upper()}
    pair2_bases = {base2.upper(), seq2[partner2].upper()}

    partner1_base = seq1[partner1]
    partner2_base = seq2[partner2]
    pair1_aln = partner1 + 1
    pair2_aln = partner2 + 1
    pair1_ung = ungapped_pos_1based(seq1, partner1)
    pair2_ung = ungapped_pos_1based(seq2, partner2)

    if pair1_bases == pair2_bases:
        typ = "No Change"
    else:
        v1, v2 = valid_pair(*pair1_bases), valid_pair(*pair2_bases)
        if v1 and v2:
            changed_site = (base1.upper() != base2.upper())
            changed_partner = (partner1_base.upper() != partner2_base.upper())
            typ = "CBC" if (changed_site and changed_partner) else "hCBC"
        elif v1 or v2:
            typ = "hCBC"
        else:
            typ = "No Pair"

    return {
        "s1_base": base1, "s1_struct": struct1[i], "s1_aln": aln_pos_1based,
        "s1_ung": ung1, "s1_pair_aln": pair1_aln, "s1_pair_ung": pair1_ung, "s1_pair_base": partner1_base,
        "s2_base": base2, "s2_struct": struct2[i], "s2_aln": aln_pos_1based,
        "s2_ung": ung2, "s2_pair_aln": pair2_aln, "s2_pair_ung": pair2_ung, "s2_pair_base": partner2_base,
        "type": typ
    }


def compact_node_info(events):
    """Extract concise 'node_X ⇒ taxon' summaries from event lines."""
    out = []
    for ev in events:
        m = re.search(r"(node_\d+).+?(Tetrix japonica|tetrix bolivari|paratettix meridionalis|tetrix bipunctata|node_\d+)", ev, flags=re.I)
        if m:
            out.append(f"{m.group(1)} ⇒ {m.group(2)}")
        else:
            out.append(ev)
    seen, uniq = set(), []
    for x in out:
        if x not in seen:
            uniq.append(x)
            seen.add(x)
    return "; ".join(uniq)


def parse_pair_arg(pairs_arg, taxa_names):
    if pairs_arg is None:
        return None
    taxa_set = set(taxa_names)
    if pairs_arg.strip().lower() == 'all':
        pairs = []
        sorted_names = sorted(taxa_names)
        for i in range(len(sorted_names)):
            for j in range(i + 1, len(sorted_names)):
                pairs.append((sorted_names[i], sorted_names[j]))
        return pairs
    # Parse comma or space separated pairs like A:B, C:D
    # Split by comma or space, then filter empty
    raw_pairs = re.split(r"[,\s]+", pairs_arg.strip())
    pairs = []
    for item in raw_pairs:
        if not item:
            continue
        if ':' not in item:
            raise ValueError(f"Invalid pair format (missing ':'): {item}")
        a, b = item.split(':', 1)
        a = a.strip()
        b = b.strip()
        if a not in taxa_set:
            raise ValueError(f"Taxon '{a}' not found in XFasta taxa.")
        if b not in taxa_set:
            raise ValueError(f"Taxon '{b}' not found in XFasta taxa.")
        pairs.append((a, b))
    return pairs

# =========================
# Main
# =========================


def main():
    p = argparse.ArgumentParser(description="Report CBC/hCBC from PAUP change list with ungapped positions.")
    p.add_argument("--changes", required=True, help="PAUP change list file")
    p.add_argument("--xfasta", required=True, help="Aligned XFasta file with sequence and structure")

    # `--pairs` is optional; we validate later that either --pairs or both --seq1/--seq2 are given
    group = p.add_mutually_exclusive_group(required=False)
    group.add_argument("--pairs", help="Pairs to compare: 'all' or comma-separated list like T1:T2,T1:T3")

    p.add_argument("--seq1", help="First taxon name in XFasta")
    p.add_argument("--seq2", help="Second taxon name in XFasta")
    p.add_argument("--out", default="cbc_results.tsv", help="Output TSV file")

    args = p.parse_args()

    # Validate that either --pairs is given or both --seq1 and --seq2 are given
    if args.pairs is None:
        if not (args.seq1 and args.seq2):
            sys.exit("Either --pairs or both --seq1 and --seq2 must be specified.")

    taxa = parse_xfasta(args.xfasta)

    if args.pairs is not None:
        pairs = parse_pair_arg(args.pairs, list(taxa.keys()))
    else:
        if args.seq1 not in taxa or args.seq2 not in taxa:
            sys.exit(f"{args.seq1} and/or {args.seq2} not found in XFasta. Found: {', '.join(taxa)}")
        pairs = [(args.seq1, args.seq2)]

    pos_events = parse_changes_with_events(args.changes)
    positions = sorted(pos_events.keys())
    if not positions:
        sys.exit("No positions found in PAUP list.")

    rows = []

    def process_pair(name1, name2):
        seq1, struct1 = taxa[name1]
        seq2, struct2 = taxa[name2]
        if len(seq1) != len(seq2):
            sys.exit(f"Alignment lengths are not equal for pair {name1} and {name2}.")

        # Precompute cached structures for performance
        mask1 = mask_structure_by_gaps(seq1, struct1)
        mask2 = mask_structure_by_gaps(seq2, struct2)
        pairs1 = find_pairs(mask1)
        pairs2 = find_pairs(mask2)
        ung1_prefix = build_ungapped_prefix(seq1)
        ung2_prefix = build_ungapped_prefix(seq2)

        for pos in positions:
            if pos < 1 or pos > len(seq1):
                continue
            res = classify_cbc_at_cached(pos, seq1, struct1, seq2, struct2,
                                         mask1, mask2, pairs1, pairs2, ung1_prefix, ung2_prefix)
            events = pos_events.get(pos, [])

            change_short = ""
            if events:
                ev0 = normalize_arrow(events[0])
                m = re.search(r"\b([ACGTU-])\s*(?:->)\s*([ACGTU-])\b", ev0)
                if m:
                    change_short = f"{normalize_base(m.group(1))}->{normalize_base(m.group(2))}"
                else:
                    change_short = ev0

            rows.append([
                name1, name2,
                pos,
                res["s1_base"], res["s1_struct"], res["s1_ung"],
                res["s1_pair_aln"], res["s1_pair_ung"], res["s1_pair_base"],
                res["s2_base"], res["s2_struct"], res["s2_ung"],
                res["s2_pair_aln"], res["s2_pair_ung"], res["s2_pair_base"],
                res["type"], change_short, "; ".join(events)
            ])

    for pair in pairs:
        process_pair(*pair)

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "seq1_name", "seq2_name",
            "aligned_pos",
            "seq1_base", "seq1_struct", "seq1_ungapped_pos", "seq1_pair_aligned_pos", "seq1_pair_ungapped_pos",
            "seq1_pair_base",
            "seq2_base", "seq2_struct", "seq2_ungapped_pos", "seq2_pair_aligned_pos", "seq2_pair_ungapped_pos",
            "seq2_pair_base",
            "type", "change", "events"
        ])
        w.writerows(rows)

    print(f"✔ Output written: {args.out} ({len(rows)} rows, {len(pairs)} pair(s))")


if __name__ == "__main__":
    main()
