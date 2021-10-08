class DFA_regular_converter {
private:
    int number_of_vertices;
    int start_vertex;
    int final_state;

    struct transition;

    mutable std::map<int, std::vector<transition> > transitions;

    void remove_double_edges();

    static std::string add(std::string a, std::string b);

    static bool available_for_mul(const std::string& s);

    static bool available_for_klini(const std::string& s);

    static std::string mul(std::string a, std::string b);

    static std::string klini(std::string a);

    void remove_vertex(int v);

    std::string get_final_regular() const;

public:

    DFA_regular_converter(const NKA& g);

    void add_edge(int v, int to, std::string s);

    std::string get_regular();
};
