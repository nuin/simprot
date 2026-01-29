# Test Data Files

This directory contains example input files for SIMPROT.

## Tree Files

Phylogenetic trees in Newick format with various numbers of taxa:

| File | Taxa | Description |
|------|------|-------------|
| `bigree40.txt` | 40 | Small test tree |
| `bigree80.txt` | 80 | Medium test tree |
| `bigree100.txt` | 100 | Medium test tree |
| `bigree111.txt` | 111 | Medium test tree |
| `bigree120.txt` | 120 | Medium test tree |
| `bigree160.txt` | 160 | Large test tree |
| `bigree320.txt` | 320 | Large test tree |

### Tree Format

Trees must be:
- **Bifurcating** (exactly two children per internal node)
- **Newick format** with branch lengths

Example:
```
((seq1:0.1,seq2:0.1):0.2,(seq3:0.1,seq4:0.1):0.2);
```

The trees use naming convention:
- Leaf nodes: `seq1`, `seq2`, etc.
- Internal nodes: assigned automatically during parsing

## Indel Distribution File

`indels.txt` - Custom indel length distribution

Format: One probability per line, representing P(length = line_number)

```
0.5     # P(length = 1)
0.25    # P(length = 2)
0.125   # P(length = 3)
...
```

Use with `-u indels.txt` option.

## Correlated Sites File

`pairs` - Correlated mutation pairs

Format: `position_i position_j correlation_value`

```
5 10 0.8
12 25 0.6
```

When position i mutates, position j has the specified probability of copying i's new residue.

Use with `-h pairs` option. Note: Indels are disabled when using correlation mode.

## Usage Examples

```bash
# Basic simulation with 40-taxon tree
simprot -f data/bigree40.txt -r 100 -s output.fasta

# Larger tree with alignment output
simprot -f data/bigree320.txt -r 200 -s seqs.fasta -a align.fasta

# Custom indel distribution
simprot -f data/bigree40.txt -r 100 -u data/indels.txt -s output.fasta

# Correlated mutations (no indels)
simprot -f data/bigree40.txt -r 100 -h data/pairs -s output.fasta
```
