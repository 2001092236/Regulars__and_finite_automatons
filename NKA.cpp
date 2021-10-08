#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <queue>
#include "constants.cpp"
#include "NKA.h"

bool NKA::dfs_deleting_epsilon_transitions(int start, int current_vertex, NKA& nka, std::vector<bool>& used) const {
    bool achieve_F = false;
    used[current_vertex] = true;

    for (auto to_vertex: transitions[{current_vertex, Constants::EPSILON}]) {
        if (!used[to_vertex]) {
            achieve_F |= dfs_deleting_epsilon_transitions(start, to_vertex, nka, used);
        }
    }

    for (char symbol: Constants::alphabet) {
        for (int to_vertex: transitions[{current_vertex, symbol}]) {
            nka.transitions[{start, symbol}].push_back(to_vertex);
        }
    }

    if (is_final_state[current_vertex]) {
        achieve_F = true;
    }
    return achieve_F;
}

void NKA::dfs_reach(int vertex, std::vector<bool>& used) const {
    used[vertex] = true;
    for (char symbol: Constants::full_alphabet) {
        for (int to_vertex: transitions[{vertex, symbol}]) {
            if (!used[to_vertex]) {
                dfs_reach(to_vertex, used);
            }
        }
    }
}

std::vector<int> NKA::get_new_class(int vertex, const std::vector<int>& prev_classes) const {
    std::vector<int> new_classes;
    new_classes.push_back(prev_classes[vertex]);
    for (char symbol: Constants::alphabet) {
        for (int to_vertex: transitions[{vertex, symbol}]) {
            new_classes.push_back(prev_classes[to_vertex]);
        }
    }
    return new_classes;
}

std::vector<std::vector<int>> NKA::get_class_matrix(const std::vector<int>& prev_classes) const {
    std::vector<std::vector<int> > class_matrix(number_of_vertices);
    for (int vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        class_matrix[vertex_id] = get_new_class(vertex_id, prev_classes);
    }
    return class_matrix;
}


void NKA::set_final_states(const std::vector<int>& new_final_states) {
    final_states = new_final_states;
    is_final_state.clear();
    is_final_state.resize(number_of_vertices, false);
    for (auto final_state: new_final_states) {
        is_final_state[final_state] = true;
    }
}

void NKA::swap(NKA& nka) {
    std::swap(number_of_vertices, nka.number_of_vertices);
    std::swap(start_vertex, nka.start_vertex);
    std::swap(final_states, nka.final_states);
    std::swap(is_final_state, nka.is_final_state);
    std::swap(transitions, nka.transitions);
}

void NKA::delete_double_edges() {
    for (int vertex_id = 0; vertex_id < number_of_vertices; ++vertex_id) {
        for (char symbol: Constants::full_alphabet) {
            std::vector<int> &to = transitions[{vertex_id, symbol}];
            sort(to.begin(), to.end());
            transitions[{vertex_id, symbol}].resize(unique(to.begin(), to.end()) - to.begin());
        }
    }
}

NKA::NKA(int number_of_vertices, int start_vertex, const std::vector<int>& final_states):
        number_of_vertices(number_of_vertices), start_vertex(start_vertex), final_states(final_states) {
    set_final_states(final_states);
}

void NKA::delete_unused_vertices() {
    std::vector<std::pair<int, int> > id_unused_vertices;

    for (auto transition: transitions) {
        if (transition.second.size() == 0) {
            id_unused_vertices.push_back(transition.first);
        }
    }

    for (auto vertex_id: id_unused_vertices) {
        transitions.erase(vertex_id);
    }
}

bool NKA::operator==(NKA& nka) {
    delete_unused_vertices();
    nka.delete_unused_vertices();
    if (number_of_vertices == nka.number_of_vertices && start_vertex == nka.start_vertex &&
                                    is_final_state == nka.is_final_state && nka.transitions == transitions) {
        return true;
    }
    return false;
}

NKA& NKA::operator=(const NKA& nka) {
    NKA copy = nka;
    swap(copy);
    return *this;
}

void NKA::add_edge(int vertex_from, int vertex_to, char symbol) {
    transitions[{vertex_from, symbol}].push_back(vertex_to);
}

NKA NKA::delete_epsilon_transitions() const {
    NKA ans = *this;

    std::vector<int> new_final_states;
    for (int vertex = 0; vertex < number_of_vertices; ++vertex) {
        std::vector<bool> used(number_of_vertices, false);
        if (dfs_deleting_epsilon_transitions(vertex, vertex, ans, used)) {
            new_final_states.push_back(vertex);
        }
    }
    ans.set_final_states(new_final_states);
    return ans;
}

NKA NKA::build_DKA() const {

    NKA dka = NKA((1LL << number_of_vertices), (1LL << start_vertex));
    std::queue<uint64_t> states;
    states.push((1LL << start_vertex));
    std::vector<bool> used((1LL << number_of_vertices), false);
    used[(1LL << start_vertex)] = true;
    while (!states.empty()) {
        uint64_t mask = states.front();
        used[mask] = true;
        states.pop();

        for (char symbol: Constants::full_alphabet) {
            uint64_t new_vertex_to = 0;
            for (int vertex = 0; vertex < number_of_vertices; ++vertex) {
                if ((mask >> vertex) & 1) {
                    for (int vertex_to: transitions[{vertex, symbol}]) {
                        new_vertex_to |= (1LL << vertex_to);
                    }
                }
            }
            if (new_vertex_to == 0) {
                continue;
            }
            dka.transitions[{mask, symbol}].push_back(new_vertex_to);
            if (!used[new_vertex_to]) {
                states.push(new_vertex_to);
            }
            used[new_vertex_to] = true;
        }
    }

    std::vector<int> new_final_states;
    uint64_t mask_of_final_states = 0;

    for (int final_state: final_states) {
        mask_of_final_states |= (1LL << final_state);
    }

    for (int new_vertex = 1; new_vertex < (1LL << number_of_vertices); ++new_vertex) {
        if (mask_of_final_states & new_vertex) {
            new_final_states.push_back(new_vertex);
        }
    }

    dka.set_final_states(new_final_states);
    return dka;
}

void NKA::transit_edges(const NKA& from, const std::vector<bool>& used, const std::vector<int>& new_numbers_of_vertices) {
    for (int vertex = 0; vertex < from.number_of_vertices; ++vertex) {
        if (!used[vertex]) {
            continue;
        }

        for (char symbol: Constants::full_alphabet) {
            for (int vertex_to: from.transitions[{vertex, symbol}]) {
                if (!used[vertex_to]) {
                    continue;
                }
                transitions[{new_numbers_of_vertices[vertex], symbol}].push_back(new_numbers_of_vertices[vertex_to]);
            }
        }
    }
}

NKA NKA::delete_unreachable_vertices() const {
    std::vector<bool> used(number_of_vertices, false);
    dfs_reach(start_vertex, used);
    std::vector<int> new_numbers_of_vertices(number_of_vertices, -1);
    new_numbers_of_vertices[start_vertex] = 0;
    int prev = 0;
    for (int vertex = 0; vertex < number_of_vertices; ++vertex) {
        if (new_numbers_of_vertices[vertex] != -1 || !used[vertex]) {
            continue;
        }
        ++prev;
        new_numbers_of_vertices[vertex] = prev;
    }

    NKA ans(prev + 1, 0);

    ans.transit_edges(*this, used, new_numbers_of_vertices);

    std::vector<int> new_final_states;

    for (int final_state: final_states) {
        if (!used[final_state]) {
            continue;
        }
        new_final_states.push_back(new_numbers_of_vertices[final_state]);
    }
    ans.set_final_states(new_final_states);
    return ans;
}

NKA NKA::get_PDKA() const {
    NKA pdka = *this;
    ++pdka.number_of_vertices;
    for (int vertex = 0; vertex <= number_of_vertices; ++vertex) {
        for (char symbol: Constants::alphabet) {
            if (pdka.transitions[{vertex, symbol}].size() == 0) {
                pdka.transitions[{vertex, symbol}].push_back(number_of_vertices);
            }
        }
    }
    return pdka;
}

NKA NKA::get_inverted_PDKA () const {
    std::vector<int> new_final_states;
    new_final_states.push_back(number_of_vertices);
    for (int vectex = 0; vectex < number_of_vertices; ++vectex) {
        if (!is_final_state[vectex]) {
            new_final_states.push_back(vectex);
        }
    }

    NKA inverted_pdka(number_of_vertices + 1, start_vertex, new_final_states);
    inverted_pdka.transitions = transitions;
    for (int vertex = 0; vertex <= number_of_vertices; ++vertex) {
        for (char symbol: Constants::alphabet) {
            if (inverted_pdka.transitions[{vertex, symbol}].size() == 0) {
                inverted_pdka.transitions[{vertex, symbol}].push_back(number_of_vertices);
            }
        }
    }
    inverted_pdka.set_final_states(new_final_states);
    return inverted_pdka;
}

std::vector<int> NKA::get_new_classes(std::vector<int> prev_classes) const {
    std::vector<int> new_classes = prev_classes;

    do{
        prev_classes = new_classes;
        std::vector<std::vector<int>> class_matrix = get_class_matrix(prev_classes);
        std::map<std::vector<int>, int> new_numbers_of_classes;
        int sz = 1;
        new_classes.clear();
        new_classes.resize(number_of_vertices, -1);
        for (int vertex = 0; vertex < number_of_vertices; ++vertex) {
            if (new_numbers_of_classes[class_matrix[vertex]] == 0) {
                new_numbers_of_classes[class_matrix[vertex]] = sz;
                new_classes[vertex] = sz - 1;
                sz++;
            } else {
                new_classes[vertex] = new_numbers_of_classes[class_matrix[vertex]] - 1;
            }
        }
    }while (new_classes != prev_classes);
    return new_classes;
}

NKA NKA::get_minimal_PDKA() const {
    std::vector<int> prev_classes(number_of_vertices, 1);
    for (int final_state: final_states) {
        prev_classes[final_state] = 0;
    }
    std::vector<int> new_classes = get_new_classes(prev_classes);

    int n_classes = *max_element(new_classes.begin(), new_classes.end()) + 1;

    std::vector<int> new_final_states;
    for (int vertex: final_states) {
        new_final_states.push_back(new_classes[vertex]);
    }
    sort(new_final_states.begin(), new_final_states.end());
    new_final_states.resize(unique(new_final_states.begin(), new_final_states.end()) - new_final_states.begin());

    NKA minimal_pdka(n_classes, new_classes[start_vertex], new_final_states);

    for (int vertex = 0; vertex < number_of_vertices; ++vertex) {
        for (char symbol: Constants::alphabet) {
            for (int vertex_to: transitions[{vertex, symbol}]) {
                minimal_pdka.add_edge(new_classes[vertex], new_classes[vertex_to], symbol);
            }
        }
    }

    minimal_pdka.delete_double_edges();
    return minimal_pdka;
}
