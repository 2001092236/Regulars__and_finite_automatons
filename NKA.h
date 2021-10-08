class NKA {
    friend class mDKA_to_regular;

public:
    bool dfs_deleting_epsilon_transitions(int start, int current_vertex, NKA& g, std::vector<bool>& used) const;

    void dfs_reach(int vertex, std::vector<bool>& used) const;

    std::vector<int> get_new_class(int vertex, const std::vector<int>& prev_classes) const;

    std::vector<std::vector<int>> get_class_matrix(const std::vector<int>& prev_classes) const;

    int number_of_vertices;
    int start_vertex;
    std::vector<int> final_states;
    std::vector<bool> is_final_state;
    mutable std::map<std::pair<int, char>, std::vector<int> > transitions;

    void set_final_states(const std::vector<int>& F1);

    void swap(NKA& g);

public:

    void delete_double_edges();

    NKA(int number_of_vertices, int start_state, const std::vector<int>& final_states = {});

    void delete_unused_vertices();

    bool operator==(NKA& g);

    NKA(const NKA& g) = default;

    NKA& operator=(const NKA& g);

    void add_edge(int v, int to, char x);

    NKA delete_epsilon_transitions() const;

    NKA build_DKA() const;

    void transit_edges(const NKA& from, const std::vector<bool>& used, const std::vector<int>& up);

    NKA delete_unreachable_vertices() const;

    NKA get_PDKA() const;

    NKA get_inverted_PDKA () const;

    std::vector<int> get_new_classes(std::vector<int> prev_classes) const;

    NKA get_minimal_PDKA() const;
};
