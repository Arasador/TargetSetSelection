class Cluster {
  int N;
  vector<vector<int>> adjacency_list;
public:
  int cluster_count;
  vector<int> clusters_separation;
  // HCS clustering
  Cluster(const vector<vector<int>>& _adjacency_list);
  vector<int> min_cut(vector<bool> vertices);
  void HCS(vector<bool> vertices);

};

Cluster::Cluster(const vector<vector<int>>& _adjacency_list) {
    adjacency_list = _adjacency_list;
    N = adjacency_list.size();
    vector<bool> initial_vertices(N, true);
    cluster_count = 0;
    clusters_separation = vector<int>(N, -1);
    HCS(initial_vertices);
}

vector<int> min_cut(const vector<int>& vertices) {
  
}

void separate_by_cut(const vector<bool>& cut,const vector<int>& vertices,
  vector<int>& H1, vector<int>& H2) {
  deque<int> queue = {vertices[0]};
  while(! queue.empty) {

  }
}

void Cluster::HCS(const vector<int>& vertices) {
  vector<int> cut = min_cut(vertices);
  int n = vertices.size();
  if (2 * cut.size() <= n) {
    for (auto v: vertices) {
      clusters[v] = cluster_count;
    }
    cluster_count ++;
    return;
  }
  vector<int> H1, H2;
  separate_by_cut(cut, vertices, H1, H2);
  HCS(H1);
  HCS(H2);
}
