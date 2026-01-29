# SIMPROT - Protein Sequence Evolution Simulator

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

SIMPROT simulates protein sequence evolution along phylogenetic trees. It generates multiple sequence alignments by simulating amino acid substitutions and insertion/deletion (indel) events based on empirically determined distributions.

## Overview

SIMPROT is widely used for benchmarking multiple sequence alignment algorithms. It provides:

- **Realistic evolution simulation** using PAM, JTT, or PMB substitution matrices
- **Empirically-based indel model** (Qian-Goldstein distribution)
- **Gamma-distributed rate heterogeneity** across sites
- **True alignment output** for benchmarking alignment accuracy

## Project Structure

```
simprot/
├── simprot-cpp/          # Modern C++20 implementation
│   ├── include/          # Header files
│   ├── src/              # Source files
│   ├── tests/            # Unit tests
│   └── CMakeLists.txt    # Build configuration
├── legacy/               # Original C/C++ implementation (v1.04)
│   ├── simprot.cpp       # Main source file
│   ├── eigen.h           # Substitution matrices
│   ├── random.c/h        # Wichmann-Hill RNG
│   └── makefile          # Build file
├── data/                 # Example tree files and test data
│   ├── bigree*.txt       # Example phylogenetic trees (40-320 taxa)
│   ├── indels.txt        # Custom indel distribution example
│   └── pairs             # Correlated site pairs example
└── README.md
```

## Implementations

### Modern C++20 Implementation (`simprot-cpp/`)

A complete rewrite using modern C++20 features with:

- Clean, modular architecture
- Full test coverage
- Exact reproducibility with the legacy version (same seed = identical output)
- No external dependencies (except standard library)

**Build:**
```bash
cd simprot-cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

**Run:**
```bash
./build/simprot -f ../data/bigree40.txt -r 100 -s sequences.fasta -a alignment.fasta
```

### Legacy Implementation (`legacy/`)

The original SIMPROT 1.04 implementation. Pre-compiled binaries are included for convenience.

**Build from source:**
```bash
cd legacy
make
```

**Dependencies:** `libpopt` for command-line parsing.

**Run:**
```bash
./simprot1.04_mac -f ../data/bigree40.txt -r 100 -s sequences.fasta -a alignment.fasta
```

## Usage

### Basic Usage

```bash
# Generate sequences along a tree
simprot -f tree.txt -r 100 -s sequences.fasta

# Generate sequences with true alignment
simprot -f tree.txt -r 100 -s sequences.fasta -a alignment.fasta

# Use JTT substitution model instead of default PMB
simprot -f tree.txt -r 100 -p 1 -s sequences.fasta
```

### Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-f, --tree` | Input tree file (Newick format, bifurcating) | required |
| `-r, --rootLength` | Root sequence length | 50 |
| `-s, --sequence` | Output sequences file (FASTA) | - |
| `-a, --alignment` | Output alignment file (FASTA with gaps) | - |
| `-p, --subModel` | Substitution model: 0=PAM, 1=JTT, 2=PMB | 2 |
| `-x, --alpha` | Gamma shape parameter (-1 for equal rates) | 1.0 |
| `-g, --indelFrequency` | Indel frequency | 0.03 |
| `-b, --branch` | Branch length scale multiplier | 1.0 |
| `-q, --InsDelRatio` | Insertion probability (vs deletion) | 0.5 |

### Input Format

**Tree file:** Newick format with branch lengths. Only bifurcating trees are supported.

```
((seq1:0.1,seq2:0.1):0.2,(seq3:0.1,seq4:0.1):0.2);
```

### Output Formats

**Sequences (`-s`):** FASTA format without gaps (leaf sequences only)
```
>seq1
MKTAYIAKQRQISFVKSH
>seq2
MKTAYIAKQRQISFVKSH
```

**Alignment (`-a`):** FASTA format with gaps showing true alignment
```
>seq1
MKT-AYIAKQRQISFVKSH
>seq2
MKTAAYIAKQRQIS-VKSH
```

## Algorithm

### Substitution Models

SIMPROT uses eigendecomposition of rate matrices for efficient probability calculation:

- **PAM** (Point Accepted Mutation) - Dayhoff et al.
- **JTT** (Jones-Taylor-Thornton) - Jones et al.
- **PMB** (Probability Matrix from Blocks) - Veerassamy et al.

### Indel Model

The default Qian-Goldstein model uses a 4-term exponential distribution based on empirical data:

```
P(length = k) = Σᵢ aᵢ * exp(-k/bᵢ)
```

Alternative models:
- **Benner/Zipfian:** `P(length = k) ∝ k^κ`
- **Custom:** Load from file

### Rate Heterogeneity

Site-specific evolutionary rates follow a gamma distribution with shape parameter α. Lower α values indicate greater rate variation.

## Citation

Please cite the following papers when using SIMPROT:

```bibtex
@article{nuin2006accuracy,
  title={The accuracy of several multiple sequence alignment methods for proteins},
  author={Nuin, Paulo AS and Wang, Zhouzhi and Tillier, Elisabeth RM},
  journal={BMC Bioinformatics},
  volume={7},
  pages={471},
  year={2006}
}

@article{pang2005simprot,
  title={SIMPROT: using an empirically determined indel distribution in simulations of protein evolution},
  author={Pang, Andy and Smith, Andrew D and Nuin, Paulo AS and Tillier, Elisabeth RM},
  journal={BMC Bioinformatics},
  volume={6},
  pages={236},
  year={2005}
}
```

## Authors

- Elisabeth Tillier
- Andy Pang
- Andrew Smith
- Robert Charlebois
- Paulo Nuin

## License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
