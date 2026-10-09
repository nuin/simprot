# SIMPROT C++20 Implementation

A modern, clean reimplementation of SIMPROT using C++20 features.

## Features

- **Modular architecture** with clear separation of concerns
- **No external dependencies** (standard library only)
- **Full test coverage** with Catch2
- **Exact reproducibility** with legacy version (same seed = identical output)
- **Type-safe interfaces** with strong typing and `std::optional`

## Building

### Requirements

- CMake 3.20+
- C++20 compatible compiler (GCC 11+, Clang 14+, MSVC 2022+)

### Build Commands

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

### Running Tests

```bash
./build/simprot_tests
```

## Architecture

```
simprot-cpp/
├── include/simprot/
│   ├── core/
│   │   ├── types.hpp           # Basic types and constants
│   │   ├── config.hpp          # SimulationConfig struct
│   │   └── random.hpp          # Wichmann-Hill RNG
│   ├── tree/
│   │   ├── tree_node.hpp       # TreeNode structure
│   │   └── tree_parser.hpp     # Newick parser
│   ├── sequence/
│   │   ├── mutable_sequence.hpp # Doubly-linked list for indels
│   │   └── alignment.hpp       # True alignment construction
│   ├── evolution/
│   │   ├── substitution_matrix.hpp # PAM/JTT/PMB matrices
│   │   ├── indel_model.hpp     # Indel length distributions
│   │   └── evolver.hpp         # Main evolution engine
│   └── io/
│       ├── newick_reader.hpp   # Tree file reader
│       └── fasta_writer.hpp    # Output writers
└── src/
    └── [corresponding .cpp files]
```

## Key Components

### WichmannHillRNG (`core/random.hpp`)

Exact port of the PAML random number generator for reproducibility:

```cpp
class WichmannHillRNG {
public:
    explicit WichmannHillRNG(int seed);
    double uniform();           // Uniform [0,1)
    double gamma(double alpha); // Gamma distribution
};
```

### SimulationConfig (`core/config.hpp`)

All simulation parameters as a structured type:

```cpp
struct SimulationConfig {
    int root_sequence_length = 50;
    SubstitutionModel substitution_model = SubstitutionModel::PMB;
    double gamma_alpha = 1.0;
    double indel_frequency = 0.03;
    double indel_ratio = 0.5;
    // ... more parameters
};
```

### TreeNode (`tree/tree_node.hpp`)

Tree structure with sequence and profile storage:

```cpp
struct TreeNode {
    std::string name;
    double distance;
    std::string sequence;
    std::vector<double> rates;
    std::unique_ptr<TreeNode> left, right;

    // Alignment profiles
    std::string left_profile_child, left_profile_self;
    std::string right_profile_child, right_profile_self;
};
```

### MutableSequence (`sequence/mutable_sequence.hpp`)

Doubly-linked list for efficient indel operations:

```cpp
class MutableSequence {
public:
    MutableSequence(const std::string& seq, const std::vector<double>& rates);

    SequenceNode* insert_before(SequenceNode* pos, size_t length, ...);
    SequenceNode* delete_at(SequenceNode* start, size_t count);

    void normalize_rates();
    std::string to_string() const;
    std::vector<double> rates() const;
};
```

### SubstitutionMatrix (`evolution/substitution_matrix.hpp`)

Interface for substitution models with eigendecomposition:

```cpp
class SubstitutionMatrix {
public:
    virtual AminoAcidIndex sample_substitution(
        AminoAcidIndex from, double t, WichmannHillRNG& rng) const = 0;
    virtual AminoAcidIndex sample_from_frequencies(
        WichmannHillRNG& rng) const = 0;
};

// Factory function
std::unique_ptr<SubstitutionMatrix> create_substitution_matrix(
    SubstitutionModel model);
```

### SequenceEvolver (`evolution/evolver.hpp`)

Main evolution engine:

```cpp
class SequenceEvolver {
public:
    SequenceEvolver(const SimulationConfig& config, WichmannHillRNG& rng);

    void init_root_sequence(TreeNode& root);
    void evolve(TreeNode& root);
};
```

## Usage Example

```cpp
#include "simprot/core/random.hpp"
#include "simprot/core/config.hpp"
#include "simprot/tree/tree_parser.hpp"
#include "simprot/evolution/evolver.hpp"
#include "simprot/io/fasta_writer.hpp"

int main() {
    // Setup
    WichmannHillRNG rng(12345);
    SimulationConfig config;
    config.root_sequence_length = 100;
    config.substitution_model = SubstitutionModel::PMB;

    // Parse tree
    NewickParser parser;
    auto root = parser.parse_file("tree.txt", rng);

    // Evolve
    SequenceEvolver evolver(config, rng);
    evolver.init_root_sequence(*root);
    evolver.evolve(*root);

    // Output
    FastaWriter writer;
    writer.write_sequences("sequences.fasta", *root);

    return 0;
}
```

## Compatibility Notes

The C++20 implementation produces **identical output** to the legacy version when given the same seed. `tests/legacy_parity/run.sh` checks this. It runs both programs over 21 combinations of trees and options (all three matrices, the three indel length models, gamma and equal rates, no indels, branch scaling, variable branches and extinction) and compares the outputs. This required careful matching of:

1. **RNG consumption order** - Same sequence of random number calls, including one draw per site when the indel frequency is 0, and the dummy "correlated mutation" draw for every profile column after the first (deleted columns count)
2. **Insertion generation** - Two-phase (all residues, then all rates)
3. **Indel position CDF** - Special handling for insertion vs deletion
4. **Zero-distance branches** - Rates not inherited, regenerated fresh
5. **Indel length CDF** - Lengths 1 .. max-1 (max = 5% of the sequence), the density cut at `DBL_EPSILON`, and the original binary search; a custom distribution always spans 200 entries and uses the unscaled branch length
6. **Insertion and deletion placement** - An insertion at position 0 is a prefix; a deletion that would run past the end is moved back to end on the last residue
7. **Amino acid frequencies** - Taken from the eigenvector row with the largest eigenvalue (row 10 for PAM, row 0 for JTT and PMB), with the cumulative table ending in 1.0
8. **Floating-point order** - `100 * distance * rate` and rate normalisation (divide by the sum, then multiply by the length) are evaluated as in the original, because PAM amplifies last-bit differences
9. **Branch extinction** - Extinct leaves and the left child of each extinct node (`FlagTree`) are left out of the output; with variable branch lengths the gamma multiplier is drawn before the extinction test

Two legacy behaviours are deliberately **not** reproduced:

- **Truncated alignments.** Legacy sometimes writes a true alignment that does not match its own sequences: after some indels at the end of a sequence it drops residues. AddressSanitizer reports a heap-buffer-overflow in `FindGapColumn`. The C++ alignment always matches the sequences, and `compare_legacy.py` reports such cases as a legacy defect.
- **Crashes.** Legacy occasionally dies before writing its output (`Killed: 9`). `run.sh` lists those seeds separately.

**PAM:** the PAM eigendecomposition in `eigen.h` / `matrix_data.hpp` has a positive eigenvalue (0.886), so PAM runs are dominated by a few amino acids (about 56% A and 26% N). The port reproduces this exactly, but the PAM matrix data needs fixing before PAM is used for science.

See `legacy/README.md` for details on the original algorithm's quirks.

## Testing

Tests use Catch2 and cover:

- RNG sequence verification against known values
- Substitution matrix probability calculations
- Indel model distributions
- Tree parsing correctness
- Legacy parity rules (`tests/test_legacy_parity.cpp`)

Run with:
```bash
./build/simprot_tests
tests/legacy_parity/run.sh build/simprot 20   # whole runs against legacy/, 20 seeds per config
```

With recent clang, Catch2 v3.4.0 trips `-Wdeprecated-literal-operator` under `-Werror`. Configure with `-DCMAKE_CXX_FLAGS=-Wno-deprecated-literal-operator` to build the tests.
