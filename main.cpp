#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <queue>
#include "constants.cpp"
using namespace std;
#include "NKA.h"
#include "DFA_regular_converter.h"


bool NKA::dfs_deleting_epsilon_transitions(int start, int current_vertex, NKA& g, vector<bool>& used) const {
    bool achieve_F = false;
    used[current_vertex] = true;

    for (auto to: transitions[{current_vertex, Constants::EPSILON}]) {
        if (!used[to]) {
            achieve_F |= dfs_deleting_epsilon_transitions(start, to, g, used);
        }
    }

    for (char c: Constants::alphabet) {
        for (int to: transitions[{current_vertex, c}]) {
            g.transitions[{start, c}].push_back(to);
        }
    }

    if (isF[current_vertex]) {
        achieve_F = true;
    }
    return achieve_F;
}

void NKA::dfs_reach(int vertex, vector<bool>& used) const {
    used[vertex] = true;
    for (char c: Constants::full_alphabet) {
        for (int to: transitions[{vertex, c}]) {
            if (!used[to]) {
                dfs_reach(to, used);
            }
        }
    }
}

vector<int> NKA::get_new_class(int vertex, const vector<int>& prev_classes) const {
    vector<int> new_classes;
    new_classes.push_back(prev_classes[vertex]);
    for (char c: Constants::alphabet) {
        for (int to: transitions[{vertex, c}]) {
            new_classes.push_back(prev_classes[to]);
        }
    }
    return new_classes;
}

vector<vector<int>> NKA::get_class_matrix(const vector<int>& prev_classes) const {
    vector<vector<int> > class_matrix(number_of_vertices);
    for (int i = 0; i < number_of_vertices; ++i) {
        class_matrix[i] = get_new_class(i, prev_classes);
    }
    return class_matrix;
}


void NKA::set_F(const vector<int>& F1) {
    F = F1;
    isF.clear();
    isF.resize(number_of_vertices, false);
    for (auto f: F1) {
        isF[f] = true;
    }
}

void NKA::swap(NKA& g) {
    ::swap(number_of_vertices, g.number_of_vertices);
    ::swap(S, g.S);
    ::swap(F, g.F);
    ::swap(isF, g.isF);
    ::swap(transitions, g.transitions);
}

void NKA::delete_double_edges() {
    for (int v = 0; v < number_of_vertices; ++v) {
        for (char c: Constants::full_alphabet) {
            vector<int> &to = transitions[{v, c}];
            sort(to.begin(), to.end());
            transitions[{v, c}].resize(unique(to.begin(), to.end()) - to.begin());
        }
    }
}

NKA::NKA(int number_of_vertices, int start_state, const vector<int>& final_states):
        number_of_vertices(number_of_vertices), S(start_state), F(final_states) {
    set_F(F);
}

void NKA::delete_unused_vertices() {
    vector<pair<int, int> > id;
    for (auto it = transitions.begin(); it != transitions.end(); ++it) {
        if ((*it).second.size() == 0) {
            id.push_back((*it).first);
        }
    }
    for (auto i: id) {
        transitions.erase(i);
    }
}

bool NKA::operator==(NKA& g) {
    delete_unused_vertices();
    g.delete_unused_vertices();
    if (number_of_vertices == g.number_of_vertices && S == g.S && isF == g.isF && g.transitions == transitions) {
        return true;
    }
    return false;
}

NKA& NKA::operator=(const NKA& g) {
    NKA cp = g;
    swap(cp);
    return *this;
}

void NKA::add_edge(int v, int to, char x) {
    transitions[{v, x}].push_back(to);
}

NKA NKA::delete_epsilon_transitions() const {
    NKA ans = *this;

    vector<int> new_final_states;
    for (int v = 0; v < number_of_vertices; ++v) {
        vector<bool> used(number_of_vertices, false);
        if (dfs_deleting_epsilon_transitions(v, v, ans, used)) {
            new_final_states.push_back(v);
        }
    }
    ans.set_F(new_final_states);
    return ans;
}

NKA NKA::build_DKA() const {

    NKA ans = NKA((1LL << number_of_vertices), (1LL << S));
    queue<uint64_t> states;
    states.push((1LL << S));
    vector<bool> used((1LL << number_of_vertices), false);
    used[(1LL << S)] = true;
    while (!states.empty()) {
        uint64_t mask = states.front();
        used[mask] = true;
        states.pop();

        for (char c: Constants::full_alphabet) {
            uint64_t to = 0;
            for (int i = 0; i < number_of_vertices; ++i) {
                if ((mask >> i) & 1) {
                    for (int j: transitions[{i, c}]) {
                        to |= (1LL << j);
                    }
                }
            }
            if (to == 0) {
                continue;
            }
            ans.transitions[{mask, c}].push_back(to);
            if (!used[to]) {
                states.push(to);
            }
            used[to] = true;
        }
    }

    vector<int> F_new;
    uint64_t mask_of_final_states = 0;

    for (int f: F) {
        mask_of_final_states |= (1LL << f);
    }

    for (int m = 1; m < (1LL << number_of_vertices); ++m) {
        if (mask_of_final_states & m) {
            F_new.push_back(m);
        }
    }

    ans.set_F(F_new);
    return ans;
}

void NKA::transit_edges(const NKA& from, const vector<bool>& used, const vector<int>& up) {
    for (int v = 0; v < from.number_of_vertices; ++v) {
        if (!used[v]) {
            continue;
        }

        for (char c: Constants::full_alphabet) {
            for (int to: from.transitions[{v, c}]) {
                if (!used[to]) {
                    continue;
                }
                transitions[{up[v], c}].push_back(up[to]);
            }
        }
    }
}

NKA NKA::delete_unreachable_vertices() const {
    vector<bool> used(number_of_vertices, false);
    dfs_reach(S, used);
    vector<int> up(number_of_vertices, -1);
    up[S] = 0;
    int prev = 0;
    for (int i = 0; i < number_of_vertices; ++i) {
        if (up[i] != -1 || !used[i]) {
            continue;
        }
        ++prev;
        up[i] = prev;
    }

    NKA ans(prev + 1, 0);

    ans.transit_edges(*this, used, up);

    vector<int> new_final_states;

    for (int f: F) {
        if (!used[f]) {
            continue;
        }
        new_final_states.push_back(up[f]);
    }
    ans.set_F(new_final_states);
    return ans;
}

NKA NKA::get_PDKA() const {
    NKA ans = *this;
    ++ans.number_of_vertices;
    for (int v = 0; v <= number_of_vertices; ++v) {
        for (char c: Constants::alphabet) {
            if (ans.transitions[{v, c}].size() == 0) {
                ans.transitions[{v, c}].push_back(number_of_vertices);
            }
        }
    }
    return ans;
}

NKA NKA::get_inverted_PDKA () const {
    vector<int> new_final_states;
    new_final_states.push_back(number_of_vertices);
    for (int i = 0; i < number_of_vertices; ++i) {
        if (!isF[i]) {
            new_final_states.push_back(i);
        }
    }

    NKA ans(number_of_vertices + 1, S, new_final_states);
    ans.transitions = transitions;
    for (int v = 0; v <= number_of_vertices; ++v) {
        for (char c: Constants::alphabet) {
            if (ans.transitions[{v, c}].size() == 0) {
                ans.transitions[{v, c}].push_back(number_of_vertices);
            }
        }
    }
    ans.set_F(new_final_states);
    return ans;
}

vector<int> NKA::get_new_classes(vector<int> prev_classes) const {
    vector<int> new_classes = prev_classes;

    do{
        prev_classes = new_classes;
        vector<vector<int>> class_matrix = get_class_matrix(prev_classes);
        map<vector<int>, int> is;
        int sz = 1;
        new_classes.clear();
        new_classes.resize(number_of_vertices, -1);
        for (int i = 0; i < number_of_vertices; ++i) {
            if (is[class_matrix[i]] == 0) {
                is[class_matrix[i]] = sz;
                new_classes[i] = sz - 1;
                sz++;
            } else {
                new_classes[i] = is[class_matrix[i]] - 1;
            }
        }
    }while (new_classes != prev_classes);
    return new_classes;
}
NKA NKA::get_minimal_PDKA() const {
    vector<int> prev_classes(number_of_vertices, 1);
    for (int f: F) {
        prev_classes[f] = 0;
    }
    vector<int> new_classes = get_new_classes(prev_classes);

    int n_classes = *max_element(new_classes.begin(), new_classes.end()) + 1;

    vector<int> new_final_states;
    for (int v: F) {
        new_final_states.push_back(new_classes[v]);
    }
    sort(new_final_states.begin(), new_final_states.end());
    new_final_states.resize(unique(new_final_states.begin(), new_final_states.end()) - new_final_states.begin());

    NKA ans(n_classes, new_classes[S], new_final_states);

    for (int v = 0; v < number_of_vertices; ++v) {
        for (char c: Constants::alphabet) {
            for (int to: transitions[{v, c}]) {
                ans.add_edge(new_classes[v], new_classes[to], c);
            }
        }
    }

    ans.delete_double_edges();
    return ans;
}


struct DFA_regular_converter::transition {
    int to;
    string regular;

    bool operator<(const transition& t) const {
        return to < t.to;
    }

    bool operator==(const transition& t) const {
        return to == t.to && regular == t.regular;
    }
};

void DFA_regular_converter::remove_double_edges() {
    for (auto it = transitions.begin(); it != transitions.end(); ++it) {
        vector<transition>& edges = (*it).second;
        vector<transition> new_e;
        sort(edges.begin(), edges.end());
        edges.resize(unique(edges.begin(), edges.end()) - edges.begin());
        if (edges.size() != 0) {
            new_e.push_back(edges[0]);
        }
        for (int i = 1; i < edges.size(); ++i) {
            if (edges[i].to == edges[i - 1].to) {
                new_e.back().regular += "+" + edges[i].regular;
            } else {
                new_e.push_back(edges[i]);
            }
        }
        edges = new_e;
    }
}

string DFA_regular_converter::add(string a, string b) {
    if (a.size() == 0) {
        return b;
    }
    if (b.size() == 0) {
        return a;
    }
    return a + " + " + b;
}

bool DFA_regular_converter::available_for_mul(const string& s) {
    if (s.size() <= 1) {
        return true;
    }

    int balance = 0;

    for (int i = 0; i < s.size(); ++i) {
        if (s[i] == '(') {
            ++balance;
            continue;
        }

        if (s[i] == ')') {
            --balance;
            continue;
        }

        if (balance == 0 && s[i] == '+') {
            return false;
        }
    }
    return true;
}

bool DFA_regular_converter::available_for_klini(const string& s) {
    if (s.size() <= 1) {
        return true;
    }
    int bal = 0;
    for (int i = 0; i < s.size(); ++i) {
        if (s[i] == '(') {
            ++bal;
        }

        if (s[i] == ')') {
            --bal;
        }
        if (bal == 0 && i + 1 != s.size()) {
            return false;
        }
    }
    return true;
}

string DFA_regular_converter::mul(string a, string b) {
    if (a.size() == 0 || a == "1") {
        return b;
    }
    if (b.size() == 0 || b == "1") {
        return a;
    }

    if (!available_for_mul(a)) {
        a = "(" + a + ")";
    }
    if (!available_for_mul(b)) {
        b = "(" + b + ")";
    }
    return a + b;
}

string DFA_regular_converter::klini(string a) {
    if (a.size() == 0) {
        return a;
    }
    if (!available_for_klini(a)) {
        a = "(" + a + ")";
    }
    return a + "*";
}

void DFA_regular_converter::remove_vertex(int v) {
    vector<transition> in;
    vector<transition> out = transitions[v];

    for (auto it = transitions.begin(); it != transitions.end(); ++it) {
        int start = (*it).first;
        if (start == v) {
            continue;
        }
        vector<transition>& to = (*it).second;
        for (int i = 0; i < to.size(); ++i) {
            string s = to[i].regular;
            int finish = to[i].to;
            if (finish == v) {
                in.push_back({start, s});
                ::swap(to[i], to.back());
                to.pop_back();
                --i;
            }
        }
    }

    string loop = "";
    for (int i = 0; i < out.size(); ++i) {
        if (out[i].to == v) {
            loop = add(loop, out[i].regular);
            ::swap(out[i], out.back());
            out.pop_back();
            --i;
        }
    }

    for (int i = 0; i < in.size(); ++i) {
        for (int j = 0; j < out.size(); ++j) {
            if (!(in[i].to != v && out[j].to != v)) {
                continue;
            }
            string s = mul(in[i].regular, mul(klini(loop), out[j].regular));
            add_edge(in[i].to, out[j].to, s);
        }
    }

    transitions[v].clear();
    transitions.erase(v);
}

string DFA_regular_converter::get_final_regular() const {
    string a, b, c, d;
    for (auto x: transitions[S]) {
        if (x.to == S) {
            a = x.regular;
        } else {
            b = x.regular;
        }
    }

    for (auto x: transitions[F]) {
        if (x.to == F) {
            d = x.regular;
        } else {
            c = x.regular;
        }
    }


    if (b == "" && F != S) {
        return "";
    }

    string result_regular = mul(klini(a), b);
    transitions.clear();
    transitions[S].push_back({F, result_regular});
    return result_regular;
}

DFA_regular_converter::DFA_regular_converter(const NKA& g): S(g.S), F(g.F[0]), n(g.number_of_vertices) {
    if (g.F.size() > 1) {
        F = n;
        for (int f: g.F) {
            transitions[f].push_back({F, "1"});
        }
        n++;
    }
    for (int v = 0; v < g.number_of_vertices; ++v) {
        for (char c: Constants::alphabet) {
            for (int to: g.transitions[{v, c}]) {
                string s;
                s += c;
                transitions[v].push_back({to, s});
            }
        }
    }
}

void DFA_regular_converter::add_edge(int v, int to, string s) {
    transitions[v].push_back({to, s});
}

string DFA_regular_converter::get_regular() {
    DFA_regular_converter r1 = *this;
    r1.remove_double_edges();
    for (int i = 0; i < n; ++i) {
        if (i == S || i == F) {
            continue;
        }
        r1.remove_vertex(i);
        r1.remove_double_edges();
    }
    return r1.get_final_regular();
}

void get_text_automata(NKA g, string name = "g_true") {
    cout << "NKA " << name << "(" << g.number_of_vertices << ", " << g.S << ", {";
    cout << g.F[0];
    for (size_t i = 1; i < g.F.size(); ++i) {
        cout << ", " << g.F[i];
    }
    cout << "});\n";

    for (auto i: g.transitions) {
        int v = i.first.first;
        char c = i.first.second;
        vector<int> to = i.second;
        for (int j = 0; j < to.size(); j++) {
            cout << name << ".add_edge(" << v << ", " << to[j] << ", " << "\'" << c << "\');\n";
        }
    }
}
