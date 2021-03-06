#include <gtest/gtest.h>
#include "NKA.h"
#include "DFA_regular_converter.h"

NKA get_AstarBstar_automata() {
    NKA g(4, 0, {0});
    g.add_edge(0, 1, 'E');
    g.add_edge(1, 3, 'a');
    g.add_edge(1, 2, 'b');
    g.add_edge(2, 3, 'a');
    g.add_edge(3, 0, 'E');
    return g;
}

TEST(AstarBstar, Deleting_epsilons) {
    NKA g = get_AstarBstar_automata();
    NKA g1 = g.delete_epsilon_transitions().delete_unreachable_vertices();
    g1.delete_double_edges();
    NKA true_g1(4, 0, {0, 3});
    true_g1.add_edge(0, 1, 'E');
    true_g1.add_edge(0, 3, 'a');
    true_g1.add_edge(0, 2, 'b');
    true_g1.add_edge(1, 3, 'a');
    true_g1.add_edge(1, 2, 'b');
    true_g1.add_edge(2, 3, 'a');
    true_g1.add_edge(3, 0, 'E');
    true_g1.add_edge(3, 3, 'a');
    true_g1.add_edge(3, 2, 'b');
    ASSERT_TRUE(g1 == true_g1);
}

TEST(AstarBstar, DKA_build) {
    NKA g = get_AstarBstar_automata();
    NKA g1 = g.delete_epsilon_transitions().delete_unreachable_vertices();

    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    NKA true_g2(4, 0, {0, 3});
    true_g2.add_edge(0, 1, 'E');
    true_g2.add_edge(0, 3, 'a');
    true_g2.add_edge(0, 2, 'b');
    true_g2.add_edge(1, 3, 'a');
    true_g2.add_edge(1, 2, 'b');
    true_g2.add_edge(2, 3, 'a');
    true_g2.add_edge(3, 0, 'E');
    true_g2.add_edge(3, 3, 'a');
    true_g2.add_edge(3, 2, 'b');

    ASSERT_TRUE(g2 == true_g2);
}

TEST(AstarBstar, PDKA_build) {
    NKA g = get_AstarBstar_automata();
    NKA g1 = g.delete_epsilon_transitions().delete_unreachable_vertices();
    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    NKA g3 = g2.get_PDKA().delete_unreachable_vertices();
    NKA true_g3(5, 0, {0, 3});
    true_g3.add_edge(0, 1, 'E');
    true_g3.add_edge(0, 3, 'a');
    true_g3.add_edge(0, 2, 'b');
    true_g3.add_edge(1, 3, 'a');
    true_g3.add_edge(1, 2, 'b');
    true_g3.add_edge(2, 3, 'a');
    true_g3.add_edge(2, 4, 'b');
    true_g3.add_edge(3, 0, 'E');
    true_g3.add_edge(3, 3, 'a');
    true_g3.add_edge(3, 2, 'b');
    true_g3.add_edge(4, 4, 'a');
    true_g3.add_edge(4, 4, 'b');

    ASSERT_TRUE(g3 == true_g3);
}

TEST(AstarBstar, minPDKA_build) {
    NKA g = get_AstarBstar_automata();
    NKA g1 = g.delete_epsilon_transitions().delete_unreachable_vertices();
    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    NKA g3 = g2.get_PDKA().delete_unreachable_vertices();
    NKA g4 = g3.get_minimal_PDKA();
    NKA true_g4(4, 0, {0});
    true_g4.add_edge(0, 0, 'a');
    true_g4.add_edge(0, 2, 'b');
    true_g4.add_edge(1, 0, 'a');
    true_g4.add_edge(1, 2, 'b');
    true_g4.add_edge(2, 0, 'a');
    true_g4.add_edge(2, 3, 'b');
    true_g4.add_edge(3, 3, 'a');
    true_g4.add_edge(3, 3, 'b');
    ASSERT_TRUE(g4 == true_g4);
}

TEST(AstarBstar, regular_build) {
    NKA g = get_AstarBstar_automata();
    NKA g1 = g.delete_epsilon_transitions().delete_unreachable_vertices();

    NKA g2 = g1.build_DKA().delete_unreachable_vertices();

    NKA g3 = g2.get_PDKA().delete_unreachable_vertices();

    NKA g4 = g3.get_minimal_PDKA();

    DFA_regular_converter r(g4);
    ASSERT_TRUE(r.get_regular() == "(a+ba)*");
}

NKA get_Epsilon_automata() {
    NKA g(1, 0, {0});
    g.add_edge(0, 0, 'E');
    return g;
}

TEST(Epsilon, Deleting_epsilons) {
    NKA g = get_Epsilon_automata();
    NKA g1 = g.delete_epsilon_transitions();
    NKA true_g1(1, 0, {0});
    true_g1.add_edge(0, 0, 'E');
    ASSERT_TRUE(true_g1 == g1);
}

TEST(Epsilon, DKA_build) {
    NKA g = get_Epsilon_automata();
    NKA g1 = g.delete_epsilon_transitions();
    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    NKA true_g2(1, 0, {0});
    true_g2.add_edge(0, 0, 'E');
    ASSERT_TRUE(true_g2 == g2);
}

TEST(Epsilon, PDKA_build) {
    NKA g = get_Epsilon_automata();
    NKA g1 = g.delete_epsilon_transitions();
    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    NKA g3 = g2.get_PDKA().delete_unreachable_vertices();
    NKA true_g3(2, 0, {0});
    true_g3.add_edge(0, 0, 'E');
    true_g3.add_edge(0, 1, 'a');
    true_g3.add_edge(0, 1, 'b');
    true_g3.add_edge(1, 1, 'a');
    true_g3.add_edge(1, 1, 'b');
    ASSERT_TRUE(true_g3 == g3);
}

TEST(Epsilon, minPDKA_build) {
    NKA g = get_Epsilon_automata();
    NKA g1 = g.delete_epsilon_transitions();
    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    NKA g3 = g2.get_PDKA().delete_unreachable_vertices();
    NKA g4 = g3.get_minimal_PDKA();
    NKA true_g4(2, 0, {0});
    true_g4.add_edge(0, 1, 'a');
    true_g4.add_edge(0, 1, 'b');
    true_g4.add_edge(1, 1, 'a');
    true_g4.add_edge(1, 1, 'b');
    ASSERT_TRUE(true_g4 == g4);
}

NKA get_A_plus_Bstar() {
    NKA g(1, 0, {0});
    g.add_edge(0, 0, 'a');
    g.add_edge(0, 0, 'b');
    return g;
}

TEST(A_plus_Bstar, Deleting_epsilons) {
    NKA g = get_A_plus_Bstar();
    NKA g1 = g.delete_epsilon_transitions();
    NKA true_g1(1, 0, {0});
    true_g1.add_edge(0, 0, 'a');
    true_g1.add_edge(0, 0, 'a');
    true_g1.add_edge(0, 0, 'b');
    true_g1.add_edge(0, 0, 'b');
    ASSERT_TRUE(g1 == true_g1);
}

TEST(A_plus_Bstar, DKA_build) {
    NKA g = get_A_plus_Bstar();
    NKA g1 = g.delete_epsilon_transitions();
    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    NKA true_g2(1, 0, {0});
    true_g2.add_edge(0, 0, 'a');
    true_g2.add_edge(0, 0, 'b');
    ASSERT_TRUE(g2 == true_g2);
}

TEST(A_plus_Bstar, Inverse_build) {
    NKA g = get_A_plus_Bstar();
    NKA g1 = g.delete_epsilon_transitions();
    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    g2 = g2.get_PDKA();
    NKA g3 = g2.get_inverted_PDKA();
    NKA true_g3(3, 0, {2, 1});
    true_g3.add_edge(0, 0, 'a');
    true_g3.add_edge(0, 0, 'b');
    true_g3.add_edge(1, 1, 'a');
    true_g3.add_edge(1, 1, 'b');
    true_g3.add_edge(2, 2, 'a');
    true_g3.add_edge(2, 2, 'b');
    ASSERT_TRUE(g3 == true_g3);
}

TEST(A_plus_Bstar, inverted_minPDKA_build) {
    NKA g = get_A_plus_Bstar();
    NKA g1 = g.delete_epsilon_transitions();
    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    g2 = g2.get_PDKA();
    NKA g3 = g2.get_inverted_PDKA();
    NKA g4 = g3.get_minimal_PDKA();
    NKA true_g4(2, 0, {1});
    true_g4.add_edge(0, 0, 'a');
    true_g4.add_edge(0, 0, 'b');
    true_g4.add_edge(1, 1, 'a');
    true_g4.add_edge(1, 1, 'b');
    ASSERT_TRUE(true_g4 == g4);
}

TEST(A_plus_Bstar, inverted_regular_build) {
    NKA g = get_A_plus_Bstar();
    NKA g1 = g.delete_epsilon_transitions();
    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    g2 = g2.get_PDKA();
    NKA g3 = g2.get_inverted_PDKA();
    NKA g4 = g3.get_minimal_PDKA();
    DFA_regular_converter r(g4);
    ASSERT_TRUE(r.get_regular() == "");
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
