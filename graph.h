//
// Created by giuseppe on 16/11/19.
//
#include <igraph/igraph.h>
#include <string>

using namespace std;

extern igraph_t generateGeometricGraph(igraph_integer_t n, igraph_real_t radius);
extern igraph_t generateBarabasiAlbertGraph(igraph_integer_t n, igraph_real_t power, igraph_integer_t m, igraph_real_t A);
extern igraph_t generateErdosRenyiGraph(igraph_integer_t n, igraph_erdos_renyi_t type, igraph_real_t param);
extern igraph_t generateRegularGraph(igraph_integer_t n, igraph_integer_t k);
extern igraph_t generateRandomGraph(int type, int n);
extern string printGraphType(int type);
