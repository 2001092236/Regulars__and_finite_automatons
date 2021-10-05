class NKA {
    friend class mDKA_to_regular;

public:
    bool dfs_deleting_epsilon_transitions(int start, int current_vertex, NKA& g, vector<bool>& used) const;

    void dfs_reach(int vertex, vector<bool>& used) const;

    vector<int> get_new_class(int vertex, const vector<int>& prev_classes) const;

    vector<vector<int>> get_class_matrix(const vector<int>& prev_classes) const;

    int number_of_vertices;
    int S;
    vector<int> F;
    vector<bool> isF;
    mutable map<pair<int, char>, vector<int> > transitions;

    void set_F(const vector<int>& F1);

    void swap(NKA& g);

public:

    void delete_double_edges();

    NKA(int number_of_vertices, int start_state, const vector<int>& final_states = {});

    void delete_unused_vertices();

    bool operator==(NKA& g);

    NKA(const NKA& g) = default;

    NKA& operator=(const NKA& g);

    void add_edge(int v, int to, char x);

    NKA delete_epsilon_transitions() const;

    NKA build_DKA() const;

    void transit_edges(const NKA& from, const vector<bool>& used, const vector<int>& up);

    NKA delete_unreachable_vertices() const;

    NKA get_PDKA() const;

    NKA get_inverted_PDKA () const;

    vector<int> get_new_classes(vector<int> prev_classes) const;

    NKA get_minimal_PDKA() const;
};
