#include "Dinic.h"

Dinic::Dinic(int _n) : MaxFlow(_n), level(_n)
{
}

Dinic::Dinic(const Dinic &obj) : MaxFlow(obj.n), level(obj.n)
{
  for(int i=0; i<obj.n; i++)
  {
    //printf("%d: ",i);
    for(int j=0; j<obj.graph[i].size(); j++)
    {
      //printf("(%d %d) ", obj.graph[i][j].v, obj.graph[i][j].cap);
      if(obj.graph[i][j].cap!=0)  //if its not a reverse arc
        addArc(i,obj.graph[i][j].v,obj.graph[i][j].cap);
    }
    //printf("\n");
  }
}


bool Dinic::buildLevelGraph(int src, int sink){
    fill(ALL(level), -1);
    while(not q.empty()) q.pop();
    level[src] = 0;
    q.push(src);
    while(not q.empty()){
      int u = q.front();
      q.pop();
      FOREACH(e, graph[u]){
        if(not e->cap or level[e->v] != -1) continue;
        level[e->v] = level[u] + 1;
        if(e->v == sink) return true;
        q.push(e->v);
      }
    }
    return false;
}

int Dinic::blockingFlow(int u, int sink, int f){
    if(u == sink or not f) return f;
    int fu = f;
    FOREACH(e, graph[u]){
      if(not e->cap or level[e->v] != level[u] + 1) continue;
      int mincap = blockingFlow(e->v, sink, min(fu, e->cap));
      if(mincap){
        graph[e->v][e->rev].cap += mincap;
        e->cap -= mincap;
        fu -= mincap;
      }
    }
    if(f == fu) level[u] = -1;
    return f - fu;
}


void Dinic::addArc(int u, int v, int cap){
    if(u == v) return;
    _Arc e(v, cap, int(graph[v].size()));
    _Arc r(u, 0, int(graph[u].size()));
    graph[u].push_back(e);
    graph[v].push_back(r);
}

bool Dinic::modifyArc(int u, int v, int cap)
{
  for(int i=0; i<graph[u].size(); i++)
  {
    if(graph[u][i].v == v && graph[u][i].cap!=0)
    {
      graph[u][i].cap = cap;
      return true;
    }
  }
  return false;
}

void Dinic::addNode()
{ 
  n++;
  vector<_Arc> new_node;
  bool aux1;
  int aux2;
  graph.push_back(new_node);
  minCutMap.push_back(aux1);
  level.push_back(aux2);
}

int Dinic::maxFlow(int src, int sink){
    flow = 0;
    while(buildLevelGraph(src, sink))
      flow += blockingFlow(src, sink, numeric_limits<int>::max());
    return flow;
}

void Dinic::runMinCut()
{
  int i;
  FORN(i, 0, n)
    if(level[i] != -1)
      minCutMap[i] = true;
    else
      minCutMap[i] = false;
}

bool Dinic::isInMinCut(int v)
{ return minCutMap[v]; }

void Dinic::print()
{
  for(int i=0; i<n; i++)
  {
    printf("%d: ",i);
    for(int j=0; j<graph[i].size(); j++)
    {
      if(graph[i][j].cap!=0)
        printf("(%d %d) ", graph[i][j].v, graph[i][j].cap);
    }
    printf("\n");
  }
}

void Dinic::create_graphviz_image(const char *filename,vector<int>& pickups, vector<int>& deliveries) {
  char command[255];
  int result;
  sprintf(command, "neato -Tgif ./graphs/graph.txt > ./graphs/%s.gif", filename);
  create_graphviz_file(pickups,deliveries);
  result = system(command);
  result = system("rm ./graphs/graph.txt");
}

void Dinic::create_graphviz_file(vector<int>& pickups, vector<int>& deliveries)
{
  const char* colors[] = {"aquamarine", "blue", "blueviolet", "burlywood", "cadetblue", "chartreuse",
    "coral", "cornflowerblue", "crimson", "darkorange", "darksalmon", "darkviolet", "gold", "goldenrod",
    "khaki", "greenyellow", "lightgrey", "palegoldenrod", "palegreen", "orangered", "tomato", "pink", "plum"};
  ofstream myfile("./graphs/graph.txt");
  if (myfile.is_open()) {
    myfile << "graph G" << "\n" << "{\noverlap=false;\n";
    myfile << "n" << 0 << "[height = 0.3, width = 0.3, label =\"" << 0 << "\", fixedsize=\"true\"];\n";
    for (int i = 0; i < pickups.size(); i++) {
      myfile << "n" << pickups[i] << "[height = 0.3, width = 0.3, label =\"" << pickups[i] << "\", style=\"filled\", fillcolor = \"" << colors[i] << "\";fixedsize=\"true\"];\n";
      myfile << "n" << deliveries[i] << "[height = 0.3, width = 0.3, shape=\"doublecircle\", label =\"" << deliveries[i] << "\", style=\"filled\", fillcolor = \"" << colors[i] << "\";fixedsize=\"true\"];\n";
    }
    myfile << "n" << 2 * pickups.size() + 1 << "[height = 0.3, width = 0.3, shape=\"doublecircle\", label =\"" << 2 * pickups.size() + 1 << "\", fixedsize=\"true\"];\n";
    for(int i=0; i<n; i++) {
      for(int j=0; j<graph[i].size(); j++) {
        if(graph[i][j].cap!=0){
          myfile << "n" << i << "-- n" << graph[i][j].v << " [dir = \"forward\", label=\"" << graph[i][j].cap << "\", color = blue];\n";
        }
      }
    }
    myfile << "}";
  }
  myfile.close();
}

bool Dinic::searchSet(int start, vector<int>& set, vector<int>& ref, bool incS[], int n)
{
  int maxValue;
  int next;
  //bool incS[50];
  bool flag, ret;
  
  set.clear();
  set.push_back(start);
  //memset(incS, false, sizeof(incS));

  incS[start] = true;
  flag = true;
  ret = false;

  while(set.size()<=n+1 && flag)
  {
    maxValue = 0;
    next = -1;
    for(int j=0; j<graph[start].size(); j++)
    {
      if(!incS[graph[start][j].v] && maxValue < graph[start][j].cap)//return edges: cap=0
      {
        maxValue = graph[start][j].cap;
        next = graph[start][j].v;
      }
    }
    if(next != -1 && next!=2*n+1)
    {
      start = next;
      set.push_back(next);
      incS[next] = true;
      flag = false;
      for(int i=0; i<ref.size(); i++)
        if(!incS[ref[i]])
          flag = true;
      if(!flag)
       ret = true;
    }
    else
      flag = false;
  }
  return ret;
}

