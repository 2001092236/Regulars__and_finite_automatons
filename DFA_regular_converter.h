class DFA_regular_converter {
private:
    int n;
    int S;
    int F;

    struct transition;


    mutable map<int, vector<transition> > transitions;

    void remove_double_edges();

    static string add(string a, string b);

    static bool available_for_mul(const string& s);

    static bool available_for_klini(const string& s);

    static string mul(string a, string b);

    static string klini(string a);

    void remove_vertex(int v);

    string get_final_regular() const;

public:

    DFA_regular_converter(const NKA& g);

    void add_edge(int v, int to, string s);

    string get_regular();
};
