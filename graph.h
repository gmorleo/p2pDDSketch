//
// Created by giuseppe on 16/11/19.
//
#include <igraph/igraph.h>
#include <string>
#include <iostream>
#include <sstream>
#include "error.h"

using namespace std;

/**
 * @brief                   This function generates a connetected random graph using the geometric model
 * @param G_graph           Where the generated graph is stored
 * @param n                 The number of vertices in the graph
 * @param radius            Graph radius
 * @return                  0 success, -10 graph generation error
 */
extern int generateGeometricGraph(igraph_t *G_graph, igraph_integer_t n, igraph_real_t radius);

/**
 * @brief                   This function generates a connected random graph using the Barabasi-Albert model
 * @param BA_graph          Where the generated graph is stored
 * @param n                 The number of vertices in the graph
 * @param power             Power of the preferential attachment. The probability that a vertex is cited is proportional to d^power+A, where d is its degree, power and A are given by arguments. In the classic preferential attachment model power=1
 * @param m                 m = number of outgoing edges generated for each vertex
 * @param A                 A = The probability that a vertex is cited is proportional to d^power+A, where d is its degree, power and A are given by argument
 * @return                  0 success, -10 graph generation error
 */
extern int generateBarabasiAlbertGraph(igraph_t *BA_graph, igraph_integer_t n, igraph_real_t power, igraph_integer_t m, igraph_real_t A);

/**
 * @brief                   This function generates a connected random graph using the Erdos-Renyi model
 * @param ER_graph          Where the generated graph is stored
 * @param n                 The number of vertices in the graph
 * @param type              IGRAPH_ERDOS_RENYI_GNM G(n,m) graph, m edges are selected uniformly randomly in a graph with n vertices\n IGRAPH_ERDOS_RENYI_GNP G(n,p) graph, every possible edge is included in the graph with probability p
 * @return                  0 success, -10 graph generation error
 */
extern int generateErdosRenyiGraph(igraph_t *ER_graph, igraph_integer_t n, igraph_erdos_renyi_t type, igraph_real_t param);

/**
 * @brief                   This function generates a connected regular random graph
 * @param R_graph           Where the generated graph is stored
 * @param n                 The number of vertices in the graph
 * @param k                 The degree of each vertex in an undirected graph. For undirected graphs, at least one of k and the number of vertices must be even.
 * @return                  0 success, -10 graph generation error
 */
extern int generateRegularGraph(igraph_t *R_graph, igraph_integer_t n, igraph_integer_t k);

/**
 * @brief                   This function generates a random graph
 * @param random_graph      Where the generated graph is stored
 * @param type              Graph distribution:\n 1 geometric\n 2 Barabasi-Albert\n 3 Erdos-Renyi\n 4 regular (clique)
 * @param n                 The number of vertices in the graph
 * @return                  0 success, -10 graph generation error
 */
extern int generateRandomGraph(igraph_t *random_graph, int type, int n);

/**
 * @brief                   This function generates a random graph
 * @param graph             Where the generated graph is stored
 * @param peers             Number of peers
 * @param graphType         Graph distribution:\n 1 geometric\n 2 Barabasi-Albert\n 3 Erdos-Renyi\n 4 regular (clique)
 * @return                  0 success, -10 graph generation error
 */
extern igraph_t* generateGraph(int peers, int graphType);

/**
 * @brief                   This function prints the name of graph distribution.
 * @param type              Graph distribution:\n 1 geometric\n 2 Barabasi-Albert\n 3 Erdos-Renyi\n 4 regular (clique)
 */
extern string printGraphProperties(igraph_t *graph);

/**
 * @brief                   This function prints the name of graph distribution.
 * @param type              Graph distribution:\n 1 geometric\n 2 Barabasi-Albert\n 3 Erdos-Renyi\n 4 regular (clique)
 */
extern string printGraphType(int type);
