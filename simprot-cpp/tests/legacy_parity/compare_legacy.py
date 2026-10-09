"""Compare a legacy SIMPROT 1.04 run with a C++20 run after normalising format.

Usage: python3 -I compare_legacy.py LEGACY_DIR CPP_DIR

Each dir holds aln.fa, seq.fa and indel.log. Legacy pads names to 20 columns,
wraps at 80 and keeps alignment columns that are gaps in every sequence; its
indel log prints the indel size before it is set and uses %f distances. So:
- seq.fa: same names in the same order, identical ungapped residues;
- aln.fa: identical after dropping all-gap columns. Legacy sometimes writes
  an alignment that does not match its own sequences (residues inserted at
  the end of a sequence are dropped). When the legacy alignment, ungapped,
  differs from the legacy sequences while the C++ one matches them, the case
  is reported as a legacy defect and not counted as a mismatch;
- indel.log: identical order of Ins/Del events. Legacy logs a distance line
  for every branch and prints an uninitialised size, C++ logs only branches
  with indels, so distances and sizes are not compared.
Prints one line per check and exits 1 on any mismatch.
"""
import sys


def fasta(path):
    recs, name = [], None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'):
            name = line[1:].strip()
            recs.append([name, ''])
        elif name is not None:
            recs[-1][1] += line.strip()
    return recs


def drop_all_gap_columns(recs):
    if not recs:
        return recs
    n = max(len(s) for _, s in recs)
    seqs = [s.ljust(n, '-') for _, s in recs]
    keep = [i for i in range(n) if any(s[i] != '-' for s in seqs)]
    return [[name, ''.join(s[i] for i in keep)] for (name, _), s in zip(recs, seqs)]


def indel_events(path):
    out = []
    for line in open(path):
        p = line.split()
        if not p:
            continue
        if p[0] == '>distance':
            out.append(('d', round(float(p[1]), 6)))
        elif p[0] in ('Ins', 'Del'):
            out.append((p[0], int(p[1])))
    return out


def main(legacy, cpp):
    ok = True

    a, b = fasta(f'{legacy}/seq.fa'), fasta(f'{cpp}/seq.fa')
    same = [n for n, _ in a] == [n for n, _ in b] and all(x[1] == y[1] for x, y in zip(a, b))
    print(f'seq.fa  {"same" if same else "DIFFERENT"} ({len(a)} vs {len(b)} sequences)')
    ok &= same

    la, ca = fasta(f'{legacy}/aln.fa'), fasta(f'{cpp}/aln.fa')
    a, b = drop_all_gap_columns(la), drop_all_gap_columns(ca)
    same = a == b
    if same:
        print(f'aln.fa  same after dropping all-gap columns ({len(a[0][1]) if a else 0} columns)')
    else:
        seqs = dict(fasta(f'{cpp}/seq.fa'))
        def consistent(aln):
            return all(s.replace('-', '') == seqs[n] for n, s in aln if n in seqs)
        if not consistent(la) and consistent(ca):
            print('aln.fa  legacy defect: the legacy alignment does not match its own sequences; '
                  'the C++ alignment does')
        else:
            print(f'aln.fa  DIFFERENT after dropping all-gap columns '
                  f'({len(a[0][1]) if a else 0} vs {len(b[0][1]) if b else 0} columns)')
            ok = False
    a = [k for k, _ in indel_events(f'{legacy}/indel.log') if k != 'd']
    b = [k for k, _ in indel_events(f'{cpp}/indel.log') if k != 'd']
    same = a == b
    print(f'indel   {"same" if same else "DIFFERENT"} Ins/Del order ({len(a)} vs {len(b)} events)')
    ok &= same
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main(*sys.argv[1:3]))
