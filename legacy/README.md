# SIMPROT Legacy Implementation (v1.04)

This directory contains the original SIMPROT implementation in C/C++.

## Files

| File | Description |
|------|-------------|
| `simprot.cpp` | Main program (~3200 lines) containing all simulation logic |
| `eigen.h` | Pre-computed eigendecomposed substitution matrices (PAM, JTT, PMB) |
| `random.c` | Wichmann-Hill random number generator (ported from PAML) |
| `random.h` | RNG header file |
| `read.c` | File I/O utilities (not currently used) |
| `makefile` | Build configuration |
| `simprot1.04_mac` | Pre-compiled binary for macOS |
| `simprot1.04_linux` | Pre-compiled binary for Linux |

## Building from Source

### Prerequisites

- C++ compiler (g++ or clang++)
- `libpopt` library for command-line parsing

**macOS (with MacPorts):**
```bash
sudo port install popt
```

**Ubuntu/Debian:**
```bash
sudo apt-get install libpopt-dev
```

### Compilation

```bash
make
```

This creates the `simprot` executable.

### Makefile Details

The makefile expects:
- popt headers in `/opt/local/include` (macOS with MacPorts)
- popt library in `/usr/local/lib`

Modify the paths if your installation differs:
```makefile
CFLAGS = -I/your/include/path
LDFLAGS = -L/your/lib/path -lpopt -lm
```

## Architecture

### Key Data Structures

**TreeNode_ (line 145):** Represents phylogenetic tree nodes
- `sequence` - Amino acid sequence (char*)
- `rate` - Per-site evolutionary rates (double*)
- `left`, `right` - Child pointers
- `distance` - Branch length to parent
- `*_profile_child`, `*_profile_self` - Alignment tracking profiles

**SequenceNode/SequenceList:** Doubly-linked list for dynamic sequences
- Enables efficient insertion/deletion operations
- Each node stores: residue, rate, mark flag

### Core Algorithm Flow

1. `main()` → Parse args, read tree, initialize matrices
2. `ReadTree()` → Parse Newick format into TreeNode structure
3. `InitRootSequence()` → Generate or load root sequence
4. `Evolve()` → Recursive tree traversal, calling `Mutate()` at each branch
5. `Mutate()` → Apply substitutions + indels based on branch length
6. `PerformIndel()` → Insert or delete sequence segments
7. `PrintFastaFormat()`/`PrintAlignment()` → Generate output

### Key Functions

| Function | Line | Description |
|----------|------|-------------|
| `GetSubstitution()` | 1020-1037 | Eigendecomposition-based substitution |
| `GetNumIndels()` | 1330-1358 | Poisson-distributed indel count |
| `InitCumulativeIndelLength()` | 1698-1740 | Qian-Goldstein indel length distribution |
| `InitCumulativeIndelLengthBenner()` | 1753-1787 | Zipfian indel length distribution |
| `InitIndel()` | 1935-1989 | Position and length selection |
| `Mutate()` | 2100-2300 | Full mutation logic |
| `Evolve()` | 2400-2600 | Recursive tree traversal |

### Global State

Substitution matrices and probability distributions are stored in global arrays:
- `matrix[20][20][4]` - Transition probability components
- `eigmat[20]` - Eigenvalues
- `SymbolCumulativeDensity[21]` - Amino acid frequencies
- `IndelCumulativeDensity` - Indel size distribution
- `Correlation[]` - Correlated mutation pairs

## Random Number Generator

The Wichmann-Hill RNG is ported from PAML for reproducibility:

```c
double rndu(void) {
    x = x * 171 % 30269;
    y = y * 172 % 30307;
    z = z * 170 % 30323;
    return fmod(x/30269.0 + y/30307.0 + z/30323.0, 1.0);
}
```

This RNG is deterministic given the same seed, enabling reproducible simulations.

## Known Quirks

1. **Rate storage:** Rates are stored per-child, not per-parent. When a parent has NULL rates, fresh gamma rates are generated for each child independently.

2. **Zero-distance branches:** When branch length is 0, only the sequence is copied to the child, not the rates. This causes fresh rate generation in subsequent mutations.

3. **Correlated mutations:** Even when no correlation file is loaded, the RNG is consumed for each position > 0 during substitutions (the correlation check always executes).

4. **Insertion CDF:** For insertions, the first residue's rate is counted twice in the cumulative distribution, and the array size is length+2 (vs length+1 for deletions).

These behaviors are preserved in the C++20 port for exact reproducibility.
