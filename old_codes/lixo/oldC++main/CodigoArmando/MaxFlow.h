#ifndef MAXFLOW_H
#define	MAXFLOW_H

#include <iostream>
#include <vector>

using namespace std;

struct _Arc{
  int v, rev;
  int cap;
  _Arc(int v_, int cap_, int rev_) : v(v_), rev(rev_), cap(cap_) {}
};

class MaxFlow
{
    protected:
      int n;
      int flow;
      vector< vector<_Arc> > graph;
      vector<bool> minCutMap;

    public:
      MaxFlow(int n_) : n(n_), graph(n_), minCutMap(n_) {}  

      virtual void addArc(int, int, int) = 0;
      virtual int maxFlow(int, int) = 0;
      virtual void runMinCut() = 0;
      virtual bool isInMinCut(int) = 0;

      virtual ~MaxFlow() {}
};

#endif	/* MAXFLOW_H */
