#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include "constants.cpp"
#include "DFA_regular_converter.h"

struct DFA_regular_converter::transition {
    int vertex_to;
    std::string regular;

    bool operator<(const transition& edge) const {
        return vertex_to < edge.vertex_to;
    }

    bool operator==(const transition& edge) const {
        return vertex_to == edge.vertex_to && regular == edge.regular;
    }
};

void DFA_regular_converter::remove_double_edges() {
    for (auto & transit: transitions) {
        std::vector<transition>& edges = transit.second;
        std::vector<transition> new_e;
        sort(edges.begin(), edges.end());
        edges.resize(unique(edges.begin(), edges.end()) - edges.begin());
        if (edges.size() != 0) {
            new_e.push_back(edges[0]);
        }
        for (int edge_id = 1; edge_id < edges.size(); ++edge_id) {
            if (edges[edge_id].vertex_to == edges[edge_id - 1].vertex_to) {
                new_e.back().regular += "+" + edges[edge_id].regular;
            } else {
                new_e.push_back(edges[edge_id]);
            }
        }
        edges = new_e;
    }
}

std::string DFA_regular_converter::add(std::string left, std::string right) {
    if (left.size() == 0) {
        return right;
    }
    if (right.size() == 0) {
        return left;
    }
    return left + " + " + right;
}

bool DFA_regular_converter::available_for_mul(const std::string& regular) {
    if (regular.size() <= 1) {
        return true;
    }

    int balance = 0;

    for (int position_in_regular = 0; position_in_regular < regular.size(); ++position_in_regular) {
        if (regular[position_in_regular] == '(') {
            ++balance;
            continue;
        }

        if (regular[position_in_regular] == ')') {
            --balance;
            continue;
        }

        if (balance == 0 && regular[position_in_regular] == '+') {
            return false;
        }
    }
    return true;
}

bool DFA_regular_converter::available_for_klini(const std::string& regular) {
    if (regular.size() <= 1) {
        return true;
    }
    int balance = 0;
    for (int position_in_regular = 0; position_in_regular < regular.size(); ++position_in_regular) {
        if (regular[position_in_regular] == '(') {
            ++balance;
        }

        if (regular[position_in_regular] == ')') {
            --balance;
        }
        if (balance == 0 && position_in_regular + 1 != regular.size()) {
            return false;
        }
    }
    return true;
}

std::string DFA_regular_converter::mul(std::string left, std::string right) {
    if (left.size() == 0 || left == "1") {
        return right;
    }
    if (right.size() == 0 || right == "1") {
        return left;
    }

    if (!available_for_mul(left)) {
        left = "(" + left + ")";
    }
    if (!available_for_mul(right)) {
        right = "(" + right + ")";
    }
    return left + right;
}

std::string DFA_regular_converter::klini(std::string regular) {
    if (regular.size() == 0) {
        return regular;
    }
    if (!available_for_klini(regular)) {
        regular = "(" + regular + ")";
    }
    return regular + "*";
}

void DFA_regular_converter::remove_vertex(int vertex) {
    std::vector<transition> in_transitions;
    std::vector<transition> out_transitions = transitions[vertex];

    for (auto & trans : transitions) {
        int start = trans.first;
        if (start == vertex) {
            continue;
        }
        std::vector<transition>& possible_transitions = trans.second;
        for (int transition_id = 0; transition_id < possible_transitions.size(); ++transition_id) {
            std::string regular = possible_transitions[transition_id].regular;
            int finish = possible_transitions[transition_id].vertex_to;
            if (finish == vertex) {
                in_transitions.push_back({start, regular});
                std::swap(possible_transitions[transition_id], possible_transitions.back());
                possible_transitions.pop_back();
                --transition_id;
            }
        }
    }

    std::string loop = "";
    for (int transition_id = 0; transition_id < out_transitions.size(); ++transition_id) {
        if (out_transitions[transition_id].vertex_to == vertex) {
            loop = add(loop, out_transitions[transition_id].regular);
            std::swap(out_transitions[transition_id], out_transitions.back());
            out_transitions.pop_back();
            --transition_id;
        }
    }

    for (auto & in_transition_id : in_transitions) {
        for (auto & out_transition_id : out_transitions) {
            if (!(in_transition_id.vertex_to != vertex && out_transition_id.vertex_to != vertex)) {
                continue;
            }
            std::string s = mul(in_transition_id.regular, mul(klini(loop), out_transition_id.regular));
            add_edge(in_transition_id.vertex_to, out_transition_id.vertex_to, s);
        }
    }

    transitions[vertex].clear();
    transitions.erase(vertex);
}

std::string DFA_regular_converter::get_final_regular() const {
    std::string loop, regular_to_finish;
    for (auto out_transition: transitions[start_vertex]) {
        if (out_transition.vertex_to == start_vertex) {
            loop = out_transition.regular;
        } else {
            regular_to_finish = out_transition.regular;
        }
    }


    if (regular_to_finish == "" && final_state != start_vertex) {
        return "";
    }

    std::string result_regular = mul(klini(loop), regular_to_finish);
    transitions.clear();
    transitions[start_vertex].push_back({final_state, result_regular});
    return result_regular;
}

DFA_regular_converter::DFA_regular_converter(const NKA& nka): start_vertex(nka.start_vertex),
                                        final_state(nka.final_states[0]), number_of_vertices(nka.number_of_vertices) {
    if (nka.final_states.size() > 1) {
        final_state = number_of_vertices;
        for (int final_state: nka.final_states) {
            transitions[final_state].push_back({final_state, "1"});
        }
        number_of_vertices++;
    }
    for (int vertex = 0; vertex < nka.number_of_vertices; ++vertex) {
        for (char symbol: Constants::alphabet) {
            for (int to: nka.transitions[{vertex, symbol}]) {
                std::string str;
                str += symbol;
                transitions[vertex].push_back({to, str});
            }
        }
    }
}

void DFA_regular_converter::add_edge(int vertex_from, int vertex_to, std::string regular) {
    transitions[vertex_from].push_back({vertex_to, regular});
}

std::string DFA_regular_converter::get_regular() {
    DFA_regular_converter converter = *this;
    converter.remove_double_edges();
    for (int vertex = 0; vertex < number_of_vertices; ++vertex) {
        if (vertex == start_vertex || vertex == final_state) {
            continue;
        }
        converter.remove_vertex(vertex);
        converter.remove_double_edges();
    }
    return converter.get_final_regular();
}
