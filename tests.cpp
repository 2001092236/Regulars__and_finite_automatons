#include "main.cpp"
#include <gtest/gtest.h>


TEST(AstarBstar, all) {
    NKA g(4, 0, {0});
    g.add_edge(0, 1, 'E');
    g.add_edge(1, 3, 'a');
    g.add_edge(1, 2, 'b');
    g.add_edge(2, 3, 'a');
    g.add_edge(3, 0, 'E');


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
    assert(g4 == true_g4);

    DFA_regular_converter r(g4);
    ASSERT_TRUE(r.get_regular() == "(a+ba)*");
}

TEST(Epsilon, all) {
    NKA g(1, 0, {0});
    g.add_edge(0, 0, 'E');

    NKA g1 = g.delete_epsilon_transitions();
    NKA true_g1(1, 0, {0});
    true_g1.add_edge(0, 0, 'E');
    ASSERT_TRUE(true_g1 == g1);


    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    NKA true_g2(1, 0, {0});
    true_g2.add_edge(0, 0, 'E');
    ASSERT_TRUE(true_g2 == g2);


    NKA g3 = g2.get_PDKA().delete_unreachable_vertices();
    NKA true_g3(2, 0, {0});
    true_g3.add_edge(0, 0, 'E');
    true_g3.add_edge(0, 1, 'a');
    true_g3.add_edge(0, 1, 'b');
    true_g3.add_edge(1, 1, 'a');
    true_g3.add_edge(1, 1, 'b');
    ASSERT_TRUE(true_g3 == g3);

    NKA g4 = g3.get_minimal_PDKA();
    NKA true_g4(2, 0, {0});
    true_g4.add_edge(0, 1, 'a');
    true_g4.add_edge(0, 1, 'b');
    true_g4.add_edge(1, 1, 'a');
    true_g4.add_edge(1, 1, 'b');
    ASSERT_TRUE(true_g4 == g4);
}

TEST(A_plus_Bstar, Inverse) {
    NKA g(1, 0, {0});
    g.add_edge(0, 0, 'a');
    g.add_edge(0, 0, 'b');
    NKA g1 = g.delete_epsilon_transitions();
    NKA true_g1(1, 0, {0});
    true_g1.add_edge(0, 0, 'a');
    true_g1.add_edge(0, 0, 'a');
    true_g1.add_edge(0, 0, 'b');
    true_g1.add_edge(0, 0, 'b');
    assert(g1 == true_g1);

    NKA g2 = g1.build_DKA().delete_unreachable_vertices();
    NKA true_g2(1, 0, {0});
    true_g2.add_edge(0, 0, 'a');
    true_g2.add_edge(0, 0, 'b');
    assert(g2 == true_g2);

    g2 = g2.get_PDKA();

    NKA g3 = g2.get_inverted_PDKA();
    NKA true_g3(3, 0, {2, 1});
    true_g3.add_edge(0, 0, 'a');
    true_g3.add_edge(0, 0, 'b');
    true_g3.add_edge(1, 1, 'a');
    true_g3.add_edge(1, 1, 'b');
    true_g3.add_edge(2, 2, 'a');
    true_g3.add_edge(2, 2, 'b');

    assert(g3 == true_g3);

    NKA g4 = g3.get_minimal_PDKA();

    NKA true_g4(2, 0, {1});
    true_g4.add_edge(0, 0, 'a');
    true_g4.add_edge(0, 0, 'b');
    true_g4.add_edge(1, 1, 'a');
    true_g4.add_edge(1, 1, 'b');
    assert(true_g4 == g4);

    DFA_regular_converter r(g4);
    assert(r.get_regular() == "");
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
