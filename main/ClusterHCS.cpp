#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <assert.h>

using namespace std;
// a structure to represent a unweighted edge in graph
struct Edge {
    int src, dest;
};

// a structure to represent a connected, undirected
// and unweighted graph as a collection of edges.
struct SimpleGraph {
    int V, E;
    vector<Edge> edges;
    vector<bool> selected;

    SimpleGraph() {}

    SimpleGraph(const vector<vector<int>>& adjacency_list){
      V = adjacency_list.size();
      E = 0;
      selected = vector<bool>(V, true);
      for (int v = 0; v < V; v ++) {
        for (auto u: adjacency_list[v]) {
          if (v < u) {
            edges.push_back({v, u});
            E ++;
          }
        }
      }
    }

    SimpleGraph extract_subgraph(vector<bool>& component) {
      SimpleGraph new_graph;
      new_graph.V = 0;
      new_graph.E = 0;
      E = 0;
      new_graph.selected = component;
      for (int v = 0; v < new_graph.selected.size(); v ++) {
        if (new_graph.selected[v]) {
          new_graph.V ++;
          selected[v] = false; // taking this vertex off
          //components[curr_component][v] = false;
        }
      }
      vector<Edge> new_edges;
      for (auto edge: edges) {
        if (new_graph.selected[edge.src] && new_graph.selected[edge.dest]) {
          new_graph.edges.push_back(edge);
          new_graph.E ++;
        }
        // in case is selected in the curr graph
        else if (selected[edge.src] && selected[edge.dest]) {
          new_edges.push_back(edge);
          E ++;
        }
      }
      V -= new_graph.V;
      edges = new_edges;
      return new_graph;
    }
};

// A structure to represent a subset for union-Find
struct subset {
    int parent;
    int rank;
};

// Function prototypes for union-Find (These functions are defined
// after kargerMinCut() )
int Find(vector<subset>&, int );
void Union(vector<subset>&, int, int);

// A very basic implementation of Karger's randomized
// algorithm for Finding the minimum cut. Please note
// that Karger's algorithm is a Monte Carlo Randomized algo
// and the cut returned by the algorithm may not be
// minimum always
int kargerMinCut(SimpleGraph& graph, vector<bool>& new_component) {
    // Get data of given graph
    int V = graph.V, E = graph.E, N = graph.selected.size();
    vector<Edge>& edge = graph.edges;
    if (edge.empty()) return 0;
    // Allocate memory for creating V subsets.
    vector<subset> subsets(N);

    // Create V subsets with single elements
    for (int v = 0; v < N; ++ v) {
      if (graph.selected[v]) {
        subsets[v].parent = v;
        subsets[v].rank = 0;
      }
    }
    // Initially there are V vertices in
    // contracted graph
    int num_vertices = V;
    // Keep contracting vertices until there are
    // 2 vertices.
    while (num_vertices > 2) {
       // Pick a random edge
       int i = rand() % E;
       // Find num_vertices (or sets) of two corners
       // of current edge
       assert(graph.selected[edge[i].src] && graph.selected[edge[i].dest]);
       int subset1 = Find(subsets, edge[i].src);
       int subset2 = Find(subsets, edge[i].dest);
       //cout << " first sub  " << subset1 << " second sub " << subset2 << endl;
       // If two corners belong to same subset,
       // then no point considering this edge
       if (subset1 == subset2)
         continue;
       // Else contract the edge (or combine the
       // corners of edge into one vertex)
       else {
          //printf("Contracting edge %d-%d\n",
                // edge[i].src, edge[i].dest);
          num_vertices --;
          Union(subsets, subset1, subset2);
       }
    }

    // Now we have two vertices (or subsets) left in
    // the contracted graph, so count the edges between
    // two components and return the count.
    int cutedges = 0;
    for (int i = 0; i < E; i++) {
      int subset1 = Find(subsets, edge[i].src);
      int subset2 = Find(subsets, edge[i].dest);
      if (subset1 != subset2) cutedges ++;
    }
    new_component = vector<bool>(N, false);
    int subset_compar;
    bool first = true;
    for (int v = 0; v < N; v ++) {
      if (graph.selected[v]) {
        //cout << v<<" subset " << Find(subsets, v) << endl;
        if (first) {
          subset_compar = Find(subsets, v);
          first = false;
        } else {
          if (subset_compar != Find(subsets, v)) {
            new_component[v] = true;
          }
        }
      }
    }
    return cutedges;
}

// A utility function to Find set of an element i
// (uses path compression technique)
int Find(vector<subset>& subsets, int i) {
    // Find root and make root as parent of i
    // (path compression)
    if (subsets[i].parent != i)
      subsets[i].parent = Find(subsets, subsets[i].parent);

    return subsets[i].parent;
}

// A function that does union of two sets of x and y
// (uses union by rank)
void Union(vector<subset>& subsets, int x, int y) {
    int xroot = Find(subsets, x);
    int yroot = Find(subsets, y);

    // Attach smaller rank tree under root of high
    // rank tree (Union by Rank)
    if (subsets[xroot].rank < subsets[yroot].rank)
        subsets[xroot].parent = yroot;
    else if (subsets[xroot].rank > subsets[yroot].rank)
        subsets[yroot].parent = xroot;

    // If ranks are same, then make one as root and
    // increment its rank by one
    else {
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
    }
}

// Creates a graph with V vertices and E edges
SimpleGraph createGraph(int V, int E) {
    SimpleGraph graph = SimpleGraph();
    graph.V = V;
    graph.E = E;
    graph.edges = vector<Edge>(E);
    graph.selected = vector<bool>(V, true);
    return graph;
}

// cluster algorithm called highly connected subgraphs
void HCS_recursive (SimpleGraph& graph, vector<vector<bool>>& components) {
  if (graph.V <= 1) return;
  //cout << "Component: ";
  int N = graph.selected.size();
  /*for (int i = 0; i < N; i ++) {
    if (graph.selected[i]) //cout << i << " ";
  }
  cout << endl; //*/
  vector<bool> new_component;
  int connectivity = kargerMinCut(graph, new_component);
  //cout << "connectivity " << connectivity << endl;
  // this component is highly connected, so leave it alone
  if (2 * connectivity >= graph.V ) {
    components.push_back(graph.selected);
    return;
  }

  vector<bool> prev_selected = graph.selected;
  SimpleGraph new_graph = graph.extract_subgraph(new_component);

  for (int i = 0; i < N; i ++) {
    if (prev_selected[i] && ! graph.selected[i] && ! new_graph.selected[i]) {
      cout << "lost one vertex" << endl;
      exit(1);
    }
    if (graph.selected[i] && new_graph.selected[i]) {
      cout << "repeated in recursion " << i << endl;
      exit(1);
    }
  }

  HCS_recursive(graph, components);
  HCS_recursive(new_graph, components);
}


vector<vector<bool>> HCS(const vector<vector<int>>& adjacency_list) {
  SimpleGraph graph = SimpleGraph(adjacency_list);
  /*for (auto edge: graph.edges) {
    cout << edge.src << " - " << edge.dest << " ,  ";
  }
  cout << endl; //*/
  vector<vector<bool>> components;
  HCS_recursive(graph, components);

  #ifdef EXCLUDE_SMALL_CLUSTERS
    vector<vector<bool>> compressed_clusters;
    vector<bool> repeated(components[0].size(), false);
    int empty_clusters = 0;
    for (auto& cluster: components) {
      int counter = 0;
      for (int i = 0; i < cluster.size(); i ++){
        if (cluster[i]) {
          if (repeated[i]) {
            cout << "repeated vertex in cluster " << i << endl;
            exit(0);
          }
          repeated[i] = true;
          counter ++;
        }
      }
      if (counter == 0) empty_clusters ++;
      if (counter >= MIN_SIZE_CLUSTER)
        compressed_clusters.push_back(cluster);
    }
    //cout << "Reduced num clusters from " << components.size() << " to " <<
      //compressed_clusters.size() << endl;
      //cout << "empty " << empty_clusters << endl;
    return compressed_clusters;
  #endif
  return components;
}


// Driver program to test above functions
/*int main()
{
    /* Let us create following unweighted graph
        0------1
        | \    |
        |   \  |
        |     \|
        2------3
    int V = 4;  // Number of vertices in graph
    int E = 5;  // Number of edges in graph //
    SimpleGraph graph = createSimpleGraph(V, E);

    // add edge 0-1
    graph.edges[0].src = 0;
    graph.edges[0].dest = 1;

    // add edge 0-2
    graph.edges[1].src = 0;
    graph.edges[1].dest = 2;

    // add edge 0-3
    graph.edges[2].src = 0;
    graph.edges[2].dest = 3;

    // add edge 1-3
    graph.edges[3].src = 1;
    graph.edges[3].dest = 3;

    // add edge 2-3
    graph.edges[4].src = 2;
    graph.edges[4].dest = 3;

    // Use a different seed value for every run.
    srand(time(NULL));
    vector<bool> v;
    printf("\nCut found by Karger's randomized algo is %d\n",
           kargerMinCut(graph, v));

    return 0;
    rickson 
} //*/
