#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <queue>

using namespace std;

class NKA {
    friend class mDKA_to_regular;

public:
    bool dfs_deleting_E(int s, int v, NKA& g, vector<bool>& used) const {
        bool is = false;
        used[v] = true;

        for (auto to: delta[{v, 'E'}])
            if (!used[to])
                is |= dfs_deleting_E(s, to, g, used);

        for (char c: alphabet)
            for (int to: delta[{v, c}])
                g.delta[{s, c}].push_back(to);

        if (isF[v])
            is = true;
        return is;
    }

    void dfs_reach(int v, vector<bool>& used) const {
        used[v] = true;
        for (char c: full_alphabet) {
            for (int to: delta[{v, c}])
                if (!used[to])
                    dfs_reach(to, used);
        }
    }

    vector<int> get_new_class(int v, const vector<int>& prev_classes) const {
        vector<int> new_cl;
        new_cl.push_back(prev_classes[v]);
        for (char c: alphabet) {
            for (int to: delta[{v, c}]) {
                new_cl.push_back(prev_classes[to]);
            }
        }
        return new_cl;
    }

    vector<vector<int>> get_class_matrix(const vector<int>& prev_classes) const {
        vector<vector<int> > A(n);
        for (int i = 0; i < n; ++i) {
            A[i] = get_new_class(i, prev_classes);
        }
        return A;
    }

    string alphabet = "abcdefghiklmnopqrstuvwxyz";
    string full_alphabet = "abcdefghiklmnopqrstuvwxyzE";
    int n;
    int S;
    vector<int> F;
    vector<bool> isF;
    mutable map<pair<int, char>, vector<int> > delta;

    void set_F(const vector<int>& F1) {
        F = F1;
        isF.clear();
        isF.resize(n, false);
        for (auto f: F1)
            isF[f] = true;
    }

    void swap(NKA& g) {
        ::swap(n, g.n);
        ::swap(S, g.S);
        ::swap(F, g.F);
        ::swap(isF, g.isF);
        ::swap(delta, g.delta);

        ::swap(alphabet, g.alphabet);
        ::swap(full_alphabet, g.full_alphabet);
    }

public:

    void delete_double_edges() {
        for (int v = 0; v < n; ++v) {
            for (char c: full_alphabet) {
                vector<int> &to = delta[{v, c}];
                sort(to.begin(), to.end());
                delta[{v, c}].resize(unique(to.begin(), to.end()) - to.begin());
            }
        }
    }

    void print() const {
        cout << "n = " << n << " S = " << S << "\n";
        for (auto i = delta.begin(); i != delta.end(); ++i) {
            int v = (*i).first.first;
            vector<int> to = (*i).second;
            if (to.size() == 0)
                continue;
            char c = (*i).first.second;
            cout << v << " -> " << c;
            cout << " -> ";
            for (int j: to)
                cout << j;
            cout << "\n";
        }
        cout << "F: ";
        for (int j: F)
            cout << j << " ";
        cout << "\n";
    }
    NKA(int n, int S, const vector<int>& F = {}, const string& alf = "ab"): n(n), S(S), F(F), alphabet(alf) {
        full_alphabet = alphabet;
        full_alphabet += 'E';
        set_F(F);
    }

    void delete_empty() {
        vector<pair<int, int> > id;
        for (auto it = delta.begin(); it != delta.end(); ++it) {
            if ((*it).second.size() == 0)
                id.push_back((*it).first);
        }
        for (auto i: id)
            delta.erase(i);
    }

    bool operator==(NKA& g) {
        delete_empty();
        g.delete_empty();
        if (n == g.n && S == g.S && isF == g.isF && g.delta == delta) {
            return true;
        }
        return false;
    }
    NKA(const NKA& g) = default;

    NKA& operator=(const NKA& g) {
        NKA cp = g;
        swap(cp);
        return *this;
    }

    void add_edge(int v, int to, char x) {
        delta[{v, x}].push_back(to);
    }

    NKA delete_E() const {
        NKA ans = *this;

        vector<int> newF;
        for (int v = 0; v < n; ++v) {
            vector<bool> used(n, false);
            if (dfs_deleting_E(v, v, ans, used))
                newF.push_back(v);
        }
        ans.set_F(newF);
        /// ans.delete_double_edges();
        return ans;
    }

    NKA build_DKA() const { /// Е-ребра уже должны быть удалены

        ///множества вершин = маски
        NKA ans = NKA((1LL << n), (1LL << S));
        queue<uint64_t> q;
        q.push((1LL << S));
        vector<bool> used((1LL << n), false);
        used[(1LL << S)] = true;
        while (!q.empty()) {
            uint64_t mask = q.front();
            used[mask] = true;
            q.pop();

            for (char c: full_alphabet) {
                uint64_t to = 0;
                for (int i = 0; i < n; ++i) {
                    if ((mask >> i) & 1) {
                        for (int j: delta[{i, c}]) {
                            to |= (1LL << j);
                        }
                    }
                }
                if (to == 0)
                    continue;
                ans.delta[{mask, c}].push_back(to);
                if (!used[to])
                    q.push(to);
                used[to] = true;
            }
        }

        vector<int> F_new;
        uint64_t my_mask = 0;

        for (int f: F)
            my_mask |= (1LL << f);

        for (int m = 1; m < (1LL << n); ++m) {
            if (my_mask & m)
                F_new.push_back(m);
        }

        ans.set_F(F_new);
        return ans;
    }

    NKA delete_unreachable_verticles() const {
        vector<bool> used(n, false);
        dfs_reach(S, used);
        vector<int> up(n, -1);
        up[S] = 0;
        int prev = 0;
        for (int i = 0; i < n; ++i) {
            if (up[i] != -1 || !used[i])
                continue;
            up[i] = ++prev;
        }

        NKA ans(prev + 1, 0);

        ///перенести ребра

        for (int v = 0; v < n; ++v) {
            if (!used[v])
                continue;

            for (char c: full_alphabet) {
                for (int to: delta[{v, c}]) {
                    if (!used[to])
                        continue;
                    ans.delta[{up[v], c}].push_back(up[to]);
                }
            }
        }

        ///конечные состояния
        vector<int> newF;

        for (int f: F) {
            if (!used[f])
                continue;
            newF.push_back(up[f]);
        }
        ans.set_F(newF);
        return ans;
    }

    NKA get_PDKA() const { ///Уже - ДКА
        NKA ans = *this;
        ++ans.n;
        for (int v = 0; v <= n; ++v) {
            for (char c: alphabet) {
                if (ans.delta[{v, c}].size() == 0) {
                    ans.delta[{v, c}].push_back(n);
                }
            }
        }
        return ans;
    }

    NKA get_inverted_PDKA () const { ///дали ПДКА
        vector<int> Fnew;
        Fnew.push_back(n);
        for (int i = 0; i < n; ++i) {
            if (!isF[i])
                Fnew.push_back(i);
        }

        NKA ans(n + 1, S, Fnew, alphabet);
        ans.delta = delta;
        for (int v = 0; v <= n; ++v) {
            for (char c: alphabet) {
                if (ans.delta[{v, c}].size() == 0) {
                    ans.delta[{v, c}].push_back(n);
                }
            }
        }
        ans.set_F(Fnew);
        return ans;
    }

    NKA get_minimal_PDKA() const { ///дали ПДКА
        vector<int> prev_classes(n, 1);
        for (int f: F)
            prev_classes[f] = 0;

        vector<int> new_classes = prev_classes;

        do{
            prev_classes = new_classes;
            vector<vector<int>> class_matrix = get_class_matrix(prev_classes);
            map<vector<int>, int> is;
            int sz = 1;
            new_classes.clear();
            new_classes.resize(n, -1);
            for (int i = 0; i < n; ++i) {
                if (is[class_matrix[i]] == 0) {
                    is[class_matrix[i]] = sz;
                    new_classes[i] = sz - 1;
                    sz++;
                } else {
                    new_classes[i] = is[class_matrix[i]] - 1;
                }
            }
        }while (new_classes != prev_classes);

        int n_classes = *max_element(new_classes.begin(), new_classes.end()) + 1;

        vector<int> Fnew;
        for (int v: F) {
            Fnew.push_back(new_classes[v]);
        }
        sort(Fnew.begin(), Fnew.end());
        Fnew.resize(unique(Fnew.begin(), Fnew.end()) - Fnew.begin());

        NKA ans(n_classes, new_classes[S], Fnew);

        for (int v = 0; v < n; ++v) {
            for (char c: alphabet) {
                for (int to: delta[{v, c}]) {
                    ans.add_edge(new_classes[v], new_classes[to], c);
                }
            }
        }

        ans.delete_double_edges();
        return ans;
    }

};

class mDKA_to_regular {
private:
    int n;
    int S;
    int F;
    mutable map<int, vector<pair<int, string> > > delta;
    string alphabet = "abcdefghiklmnopqrstuvwxyz";
    string full_alphabet = "abcdefghiklmnopqrstuvwxyzE";

    void remove_double_edges() {
        for (auto it = delta.begin(); it != delta.end(); ++it) {
            vector<pair<int, string> >& edges = (*it).second;
            vector<pair<int, string> > new_e;
            sort(edges.begin(), edges.end());
            edges.resize(unique(edges.begin(), edges.end()) - edges.begin());
            if (edges.size() != 0)
                new_e.push_back(edges[0]);
            for (int i = 1; i < edges.size(); ++i) {
                if (edges[i].first == edges[i - 1].first) {
                    new_e.back().second += "+" + edges[i].second;
                } else {
                    new_e.push_back(edges[i]);
                }
            }
            edges = new_e;
        }
    }

    static string add(string a, string b) {
        if (a.size() == 0)
            return b;
        if (b.size() == 0)
            return a;

        return a + " + " + b;
    }

    static bool free_to_mul(const string& s) {
        if (s.size() == 0)
            return true;
        if (s.size() == 1)
            return true;
        int bal = 0;

        for (int i = 0; i < s.size(); ++i) {
            if (s[i] == '(') {
                ++bal;
                continue;
            }

            if (s[i] == ')') {
                --bal;
                continue;
            }

            if (bal == 0 && s[i] == '+')
                return false;
        }
        return true;
    }

    static bool free_to_klini(const string& s) {
        if (s.size() == 0)
            return true;
        if (s.size() == 1)
            return true;
        int bal = 0;
        for (int i = 0; i < s.size(); ++i) {
            if (s[i] == '(') {
                ++bal;
            }

            if (s[i] == ')') {
                --bal;
            }
            if (bal == 0 && i + 1 != s.size())
                return false;
        }
        return true;
    }

    static string mul(string a, string b) {
        if (a.size() == 0 || a == "1")
            return b;
        if (b.size() == 0 || b == "1")
            return a;

        if (!free_to_mul(a))
            a = "(" + a + ")";
        if (!free_to_mul(b))
            b = "(" + b + ")";
        return a + b;
    }

    static string klini(string a) {
        if (a.size() == 0)
            return a;
        if (!free_to_klini(a))
            a = "(" + a + ")";
        return a + "*";
    }

    void remove_vertex(int v) { /// v - не первая и не последняя
        vector<pair<int, string> > in;
        vector<pair<int, string> > out = delta[v];

        for (auto it = delta.begin(); it != delta.end(); ++it) {
            int start = (*it).first;
            if (start == v)
                continue;
            vector<pair<int, string> >& to = (*it).second;
            for (int i = 0; i < to.size(); ++i) {
                string s = to[i].second;
                int finish = to[i].first;
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
            if (out[i].first == v) {
                loop = add(loop, out[i].second);
                ::swap(out[i], out.back());
                out.pop_back();
                --i;
            }
        }

        for (int i = 0; i < in.size(); ++i) {
            for (int j = 0; j < out.size(); ++j) {
                if (!(in[i].first != v && out[j].first != v))
                    continue;
                string s = mul(in[i].second, mul(klini(loop), out[j].second));
                add_edge(in[i].first, out[j].first, s);
            }
        }

        delta[v].clear();
        delta.erase(v);
    }

    string get_final_regular() const {
        string a, b, c, d;
        for (auto x: delta[S]) {
            if (x.first == S) {
                a = x.second;
            } else {
                b = x.second;
            }
        }

        for (auto x: delta[F]) {
            if (x.first == F) {
                d = x.second;
            } else {
                c = x.second;
            }
        }

        ///cout << "       ERTERTERTERT:  a=" << a << "$ b=" << b << "$ c=" << c << "$ d=" << d << "$\n\n";

        if (b == "" && F != S)
            return "";


        /// string str = mul(klini(a), mul(b, klini(add(d, mul(c, mul(klini(a), b)))))); /// Из конца не могут выходить ребра
        string str = mul(klini(a), b);
        delta.clear();
        delta[S].push_back({F, str});
        return str;
    }

public:

    void print() {
        cout << "n = " << n << " S = " << S << "\n";
        for (auto i = delta.begin(); i != delta.end(); ++i) {
            int v = (*i).first;

            for (auto x: delta[v]) {
                cout << v << " -> " << x.second << " -> " << x.first << "\n";
            }
        }
        cout << "F: ";
        cout << F;
        cout << "\n";
    }
    mDKA_to_regular(const NKA& g): S(g.S), F(g.F[0]), n(g.n), alphabet(g.alphabet), full_alphabet(g.full_alphabet) {
        if (g.F.size() > 1) {
            F = n;
            for (int f: g.F)
                delta[f].push_back({F, "1"});
            n++;
        }
        for (int v = 0; v < g.n; ++v) {
            for (char c: alphabet) {
                for (int to: g.delta[{v, c}]) {
                    string s;
                    s += c;
                    delta[v].push_back({to, s});
                }
            }
        }
    }

    void add_edge(int v, int to, string s) {
        delta[v].push_back({to, s});
    }

    string get_regular() {
        mDKA_to_regular r1 = *this;
        r1.remove_double_edges();
        for (int i = 0; i < n; ++i) {
            if (i == S || i == F)
                continue;
            r1.remove_vertex(i);
            r1.remove_double_edges();
        }
        return r1.get_final_regular();
    }

};

void get_text_automata(NKA g, string name = "g_true") {
    cout << "NKA " << name << "(" << g.n << ", " << g.S << ", {";
    cout << g.F[0];
    for (size_t i = 1; i < g.F.size(); ++i)
        cout << ", " << g.F[i];
    cout << "});\n";

    for (auto i: g.delta) {
        int v = i.first.first;
        char c = i.first.second;
        vector<int> to = i.second;
        for (int j = 0; j < to.size(); j++) {
            cout << name << ".add_edge(" << v << ", " << to[j] << ", " << "\'" << c << "\');\n";
        }
    }
}