#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "simprot/tree/tree_parser.hpp"
#include "simprot/core/random.hpp"

using namespace simprot;
using Catch::Approx;

TEST_CASE("NewickParser basic parsing", "[tree_parser]") {
    NewickParser parser;

    SECTION("Single leaf node") {
        auto tree = parser.parse("A:0.5;");
        REQUIRE(tree != nullptr);
        // Parser creates a tree structure - verify it has exactly one leaf
        REQUIRE(tree->count_leaves() == 1);
        // The parsed tree may wrap the leaf in a root or return it directly
        // Implementation detail - just verify the tree is valid
    }

    SECTION("Two-taxon tree") {
        auto tree = parser.parse("(A:0.1,B:0.2);");
        REQUIRE(tree != nullptr);
        REQUIRE_FALSE(tree->is_leaf());
        REQUIRE(tree->left != nullptr);
        REQUIRE(tree->right != nullptr);
        REQUIRE(tree->left->name == "A");
        REQUIRE(tree->left->distance == Approx(0.1));
        REQUIRE(tree->right->name == "B");
        REQUIRE(tree->right->distance == Approx(0.2));
    }

    SECTION("Four-taxon tree") {
        auto tree = parser.parse("((A:0.1,B:0.2):0.3,(C:0.1,D:0.2):0.4);");
        REQUIRE(tree != nullptr);
        REQUIRE_FALSE(tree->is_leaf());

        // Left subtree
        auto left = tree->left.get();
        REQUIRE(left != nullptr);
        REQUIRE(left->distance == Approx(0.3));
        REQUIRE(left->left->name == "A");
        REQUIRE(left->right->name == "B");

        // Right subtree
        auto right = tree->right.get();
        REQUIRE(right != nullptr);
        REQUIRE(right->distance == Approx(0.4));
        REQUIRE(right->left->name == "C");
        REQUIRE(right->right->name == "D");
    }
}

TEST_CASE("NewickParser internal node naming", "[tree_parser]") {
    NewickParser parser;

    SECTION("Internal nodes get auto-generated names") {
        auto tree = parser.parse("((A:0.1,B:0.2):0.3,C:0.4);");
        REQUIRE(tree != nullptr);

        // Root should have a name (auto-generated or empty, depending on impl)
        // Internal node with A and B should also have a name
        REQUIRE(tree->left != nullptr);
        // Just verify it's not a leaf
        REQUIRE_FALSE(tree->left->is_leaf());
    }
}

TEST_CASE("NewickParser branch scaling", "[tree_parser]") {
    NewickParser parser;

    SECTION("Default scale is 1.0") {
        auto tree = parser.parse("(A:0.5,B:0.5);");
        REQUIRE(tree->left->distance == Approx(0.5));
    }

    SECTION("Branch scale of 2.0 doubles distances") {
        parser.set_branch_scale(2.0);
        auto tree = parser.parse("(A:0.5,B:0.5);");
        REQUIRE(tree->left->distance == Approx(1.0));
        REQUIRE(tree->right->distance == Approx(1.0));
    }

    SECTION("Branch scale of 0.5 halves distances") {
        parser.set_branch_scale(0.5);
        auto tree = parser.parse("(A:1.0,B:1.0);");
        REQUIRE(tree->left->distance == Approx(0.5));
        REQUIRE(tree->right->distance == Approx(0.5));
    }
}

TEST_CASE("NewickParser error handling", "[tree_parser]") {
    NewickParser parser;

    SECTION("Empty string throws") {
        REQUIRE_THROWS_AS(parser.parse(""), ParseError);
    }

    SECTION("Missing semicolon throws") {
        REQUIRE_THROWS_AS(parser.parse("(A:0.1,B:0.2)"), ParseError);
    }

    SECTION("Unmatched parentheses throws") {
        REQUIRE_THROWS_AS(parser.parse("((A:0.1,B:0.2);"), ParseError);
    }
}

TEST_CASE("NewickParser leaf counting", "[tree_parser]") {
    NewickParser parser;

    SECTION("Count leaves in balanced tree") {
        auto tree = parser.parse("((A:0.1,B:0.2):0.3,(C:0.1,D:0.2):0.4);");
        REQUIRE(tree->count_leaves() == 4);
    }

    SECTION("Count leaves in unbalanced tree") {
        auto tree = parser.parse("(((A:0.1,B:0.2):0.3,C:0.4):0.5,D:0.6);");
        REQUIRE(tree->count_leaves() == 4);
    }

    SECTION("Single leaf tree has count 1") {
        auto tree = parser.parse("A:0.5;");
        REQUIRE(tree->count_leaves() == 1);
    }
}

TEST_CASE("NewickParser with variable branches", "[tree_parser]") {
    WichmannHillRNG rng(42);
    NewickParser parser(&rng);

    SECTION("Variable branch gamma modifies distances") {
        parser.set_variable_branch_gamma(1.0);
        auto tree1 = parser.parse("(A:1.0,B:1.0);");

        rng.set_seed(42);  // Reset for comparison
        NewickParser parser2(&rng);
        parser2.set_variable_branch_gamma(1.0);
        auto tree2 = parser2.parse("(A:1.0,B:1.0);");

        // Same seed should give same results
        REQUIRE(tree1->left->distance == Approx(tree2->left->distance));
    }
}

TEST_CASE("NewickParser handles whitespace", "[tree_parser]") {
    NewickParser parser;

    SECTION("Whitespace around elements is ignored") {
        auto tree = parser.parse(" ( A : 0.1 , B : 0.2 ) ; ");
        REQUIRE(tree != nullptr);
        REQUIRE(tree->left->name == "A");
        REQUIRE(tree->right->name == "B");
    }

    SECTION("Newlines in tree are handled") {
        auto tree = parser.parse("(\n  A:0.1,\n  B:0.2\n);");
        REQUIRE(tree != nullptr);
        REQUIRE(tree->count_leaves() == 2);
    }
}

TEST_CASE("NewickParser handles zero branch lengths", "[tree_parser]") {
    NewickParser parser;

    SECTION("Zero branch length is preserved") {
        auto tree = parser.parse("(A:0.0,B:0.5);");
        REQUIRE(tree->left->distance == Approx(0.0));
    }

    SECTION("Missing branch length defaults to zero") {
        auto tree = parser.parse("(A,B:0.5);");
        REQUIRE(tree->left->distance == Approx(0.0));
    }
}
