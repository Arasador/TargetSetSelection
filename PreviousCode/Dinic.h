#ifndef DINIC_H
#define DINIC_H

#include <vector>
#include <queue>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>

#include "MaxFlow.h"

using namespace std;

//######################################
//#####  Dinic MaxFlow Algorithm ######
//######################################
//#define NDEBUG
#if defined(NDEBUG)
#define DBG_CODE(cb...)
#else
#define DBG_CODE(cb...) cb
#endif
#define WRITE(x) DBG_CODE(cout << x << endl)
#define WATCH(x) DBG_CODE(cout << #x << "=" << x << endl)
#define FORN(i, a, b) for(typeof(b) i = (a); i < (b); i++)
#define ALL(x) x.begin(), x.end()
#define FOREACH(i, c) for(typeof((c).begin()) i = (c).begin(); i != (c).end(); i++)


class Dinic : public MaxFlow
{
    private:
      vector<int> level;
      queue<int> q;
      bool buildLevelGraph(int, int);
      int blockingFlow(int, int, int);

    public:
      Dinic(int);
      Dinic(const Dinic &obj);

      virtual void addArc(int, int, int);
      virtual int maxFlow(int, int);
      virtual void runMinCut();
      virtual bool isInMinCut(int);
      
      bool modifyArc(int, int, int);
      void addNode();
      void print();
      bool searchSet(int, vector<int>&, vector<int>&, bool[], int);
      void create_graphviz_file(vector<int>&,vector<int>&);
      void create_graphviz_image(const char *,vector<int>& pickups, vector<int>& deliveries);


};

#endif  /* DINIC_H */
