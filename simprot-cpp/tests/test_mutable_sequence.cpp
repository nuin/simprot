#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "simprot/sequence/mutable_sequence.hpp"
#include "simprot/evolution/substitution_matrix.hpp"
#include "simprot/core/random.hpp"

using namespace simprot;
using Catch::Approx;

TEST_CASE("MutableSequence construction", "[mutable_sequence]") {
    SECTION("Empty sequence") {
        MutableSequence seq;
        REQUIRE(seq.empty());
        REQUIRE(seq.size() == 0);
        REQUIRE(seq.head() == nullptr);
        REQUIRE(seq.tail() == nullptr);
    }

    SECTION("From string with rates") {
        std::vector<double> rates = {1.0, 1.5, 0.5, 2.0, 1.0};
        MutableSequence seq("ACDEF", rates);

        REQUIRE(seq.size() == 5);
        REQUIRE_FALSE(seq.empty());
        REQUIRE(seq.to_string() == "ACDEF");
    }

    SECTION("Rates are preserved") {
        std::vector<double> rates = {1.0, 1.5, 0.5, 2.0, 1.0};
        MutableSequence seq("ACDEF", rates);

        auto extracted_rates = seq.rates();
        REQUIRE(extracted_rates.size() == 5);
        REQUIRE(extracted_rates[0] == Approx(1.0));
        REQUIRE(extracted_rates[1] == Approx(1.5));
        REQUIRE(extracted_rates[2] == Approx(0.5));
    }
}

TEST_CASE("MutableSequence to_string", "[mutable_sequence]") {
    std::vector<double> rates(10, 1.0);
    MutableSequence seq("ACDEFGHIKL", rates);

    SECTION("Returns correct sequence") {
        REQUIRE(seq.to_string() == "ACDEFGHIKL");
    }

    SECTION("Empty sequence returns empty string") {
        MutableSequence empty_seq;
        REQUIRE(empty_seq.to_string() == "");
    }
}

TEST_CASE("MutableSequence node_at", "[mutable_sequence]") {
    std::vector<double> rates = {1.0, 1.5, 0.5, 2.0, 1.0};
    MutableSequence seq("ACDEF", rates);

    SECTION("Returns correct nodes") {
        auto node0 = seq.node_at(0);
        REQUIRE(node0 != nullptr);
        REQUIRE(node0->residue == 'A');

        auto node2 = seq.node_at(2);
        REQUIRE(node2 != nullptr);
        REQUIRE(node2->residue == 'D');

        auto node4 = seq.node_at(4);
        REQUIRE(node4 != nullptr);
        REQUIRE(node4->residue == 'F');
    }

    SECTION("Out of range returns nullptr") {
        REQUIRE(seq.node_at(5) == nullptr);
        REQUIRE(seq.node_at(100) == nullptr);
    }
}

TEST_CASE("MutableSequence insertion", "[mutable_sequence]") {
    std::vector<double> rates = {1.0, 1.0, 1.0};
    MutableSequence seq("ACD", rates);
    WichmannHillRNG rng(12345);
    auto matrix = create_substitution_matrix(SubstitutionModel::PMB);

    SECTION("Insert at head") {
        auto before = seq.node_at(0);
        auto inserted = seq.insert_before(before, 2, rng, *matrix, 1.0);

        REQUIRE(seq.size() == 5);
        REQUIRE(inserted != nullptr);
        // Head should be one of the inserted nodes
        REQUIRE(seq.head() == inserted);
    }

    SECTION("Insert in middle") {
        auto before = seq.node_at(1);  // 'C'
        seq.insert_before(before, 2, rng, *matrix, 1.0);

        REQUIRE(seq.size() == 5);
        // Original 'C' should still be at some position
        std::string result = seq.to_string();
        REQUIRE(result.find('C') != std::string::npos);
    }

    SECTION("Insert at tail (before nullptr)") {
        seq.insert_before(nullptr, 2, rng, *matrix, 1.0);

        REQUIRE(seq.size() == 5);
    }

    SECTION("Inserted residues are valid amino acids") {
        seq.insert_before(seq.node_at(1), 10, rng, *matrix, 1.0);

        std::string result = seq.to_string();
        for (char c : result) {
            bool valid = false;
            for (char aa : kAminoAcidChars) {
                if (c == aa) {
                    valid = true;
                    break;
                }
            }
            REQUIRE(valid);
        }
    }
}

TEST_CASE("MutableSequence deletion", "[mutable_sequence]") {
    std::vector<double> rates = {1.0, 1.0, 1.0, 1.0, 1.0};
    MutableSequence seq("ACDEF", rates);

    SECTION("Delete from head") {
        auto next = seq.delete_at(seq.node_at(0), 2);

        REQUIRE(seq.size() == 3);
        REQUIRE(seq.to_string() == "DEF");
        REQUIRE(next != nullptr);
        REQUIRE(next->residue == 'D');
    }

    SECTION("Delete from middle") {
        seq.delete_at(seq.node_at(1), 2);  // Delete 'C' and 'D'

        REQUIRE(seq.size() == 3);
        REQUIRE(seq.to_string() == "AEF");
    }

    SECTION("Delete from tail") {
        seq.delete_at(seq.node_at(3), 2);  // Delete 'E' and 'F'

        REQUIRE(seq.size() == 3);
        REQUIRE(seq.to_string() == "ACD");
    }

    SECTION("Delete entire sequence") {
        seq.delete_at(seq.node_at(0), 5);

        REQUIRE(seq.size() == 0);
        REQUIRE(seq.empty());
        REQUIRE(seq.to_string() == "");
    }

    SECTION("Delete count exceeds remaining") {
        seq.delete_at(seq.node_at(3), 10);  // Try to delete more than available

        // Should only delete what's available
        REQUIRE(seq.size() == 3);
        REQUIRE(seq.to_string() == "ACD");
    }
}

TEST_CASE("MutableSequence marks", "[mutable_sequence]") {
    std::vector<double> rates = {1.0, 1.0, 1.0, 1.0, 1.0};
    MutableSequence seq("ACDEF", rates);

    SECTION("Clear marks sets all to zero") {
        // First set some marks
        seq.node_at(1)->mark = 1;
        seq.node_at(3)->mark = 2;

        seq.clear_marks();

        for (std::size_t i = 0; i < seq.size(); ++i) {
            REQUIRE(seq.node_at(i)->mark == 0);
        }
    }

    SECTION("for_each visits all nodes") {
        int count = 0;
        seq.for_each([&count](const SequenceNode&) {
            count++;
        });

        REQUIRE(count == 5);
    }

    SECTION("for_each can modify nodes") {
        seq.for_each([](SequenceNode& node) {
            node.mark = 42;
        });

        for (std::size_t i = 0; i < seq.size(); ++i) {
            REQUIRE(seq.node_at(i)->mark == 42);
        }
    }
}

TEST_CASE("MutableSequence rate normalization", "[mutable_sequence]") {
    std::vector<double> rates = {2.0, 4.0, 6.0, 8.0, 10.0};  // Sum = 30
    MutableSequence seq("ACDEF", rates);

    SECTION("Normalize makes sum equal to length") {
        seq.normalize_rates();

        auto normalized = seq.rates();
        double sum = 0.0;
        for (double r : normalized) {
            sum += r;
        }

        // Sum should equal sequence length (5)
        REQUIRE(sum == Approx(5.0).margin(1e-10));
    }

    SECTION("Normalize preserves relative rates") {
        seq.normalize_rates();

        auto normalized = seq.rates();
        // Original ratios: 2:4:6:8:10 = 1:2:3:4:5
        REQUIRE(normalized[1] == Approx(2.0 * normalized[0]).margin(1e-10));
        REQUIRE(normalized[2] == Approx(3.0 * normalized[0]).margin(1e-10));
    }
}

TEST_CASE("MutableSequence move semantics", "[mutable_sequence]") {
    std::vector<double> rates = {1.0, 1.0, 1.0};
    MutableSequence seq1("ACD", rates);

    SECTION("Move constructor") {
        MutableSequence seq2(std::move(seq1));

        REQUIRE(seq2.size() == 3);
        REQUIRE(seq2.to_string() == "ACD");
        REQUIRE(seq1.empty());  // NOLINT: testing moved-from state
    }

    SECTION("Move assignment") {
        MutableSequence seq2;
        seq2 = std::move(seq1);

        REQUIRE(seq2.size() == 3);
        REQUIRE(seq2.to_string() == "ACD");
        REQUIRE(seq1.empty());  // NOLINT: testing moved-from state
    }
}

TEST_CASE("MutableSequence linked list integrity", "[mutable_sequence]") {
    std::vector<double> rates = {1.0, 1.0, 1.0, 1.0, 1.0};
    MutableSequence seq("ACDEF", rates);

    SECTION("Forward traversal reaches all nodes") {
        int count = 0;
        for (auto* node = seq.head(); node != nullptr; node = node->next) {
            count++;
        }
        REQUIRE(count == 5);
    }

    SECTION("Backward traversal reaches all nodes") {
        int count = 0;
        for (auto* node = seq.tail(); node != nullptr; node = node->prev) {
            count++;
        }
        REQUIRE(count == 5);
    }

    SECTION("Head has no prev") {
        REQUIRE(seq.head()->prev == nullptr);
    }

    SECTION("Tail has no next") {
        REQUIRE(seq.tail()->next == nullptr);
    }

    SECTION("After insertion, list integrity maintained") {
        WichmannHillRNG rng(42);
        auto matrix = create_substitution_matrix(SubstitutionModel::PMB);
        seq.insert_before(seq.node_at(2), 3, rng, *matrix, 1.0);

        // Forward
        int fwd_count = 0;
        for (auto* node = seq.head(); node != nullptr; node = node->next) {
            fwd_count++;
        }

        // Backward
        int bwd_count = 0;
        for (auto* node = seq.tail(); node != nullptr; node = node->prev) {
            bwd_count++;
        }

        REQUIRE(fwd_count == 8);
        REQUIRE(bwd_count == 8);
        REQUIRE(seq.size() == 8);
    }

    SECTION("After deletion, list integrity maintained") {
        seq.delete_at(seq.node_at(1), 2);

        // Forward
        int fwd_count = 0;
        for (auto* node = seq.head(); node != nullptr; node = node->next) {
            fwd_count++;
        }

        // Backward
        int bwd_count = 0;
        for (auto* node = seq.tail(); node != nullptr; node = node->prev) {
            bwd_count++;
        }

        REQUIRE(fwd_count == 3);
        REQUIRE(bwd_count == 3);
        REQUIRE(seq.size() == 3);
    }
}
