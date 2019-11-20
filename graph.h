//
// Created by giuseppe on 16/11/19.
//
#include <igraph/igraph.h>
#include <string>

using namespace std;

/**
 * @brief                   This function generates a connetected random graph using the geometric model
 * @param n                 The number of vertices in the graph
 * @param radius            Graph radius
 * @return                  Connected random graph using the geometric model
 */
extern igraph_t generateGeometricGraph(igraph_integer_t n, igraph_real_t radius);

/**
 * @brief                   This function generates a connected random graph using the Barabasi-Albert model
 * @param n                 The number of vertices in the graph
 * @param power             Power of the preferential attachment. The probability that a vertex is cited is proportional to d^power+A, where d is its degree, power and A are given by arguments. In the classic preferential attachment model power=1
 * @param m                 m = number of outgoing edges generated for each vertex
 * @param A                 A = The probability that a vertex is cited is proportional to d^power+A, where d is its degree, power and A are given by argument
 * @return                  Connected random graph using the Barabasi-Albert model
 */
extern igraph_t generateBarabasiAlbertGraph(igraph_integer_t n, igraph_real_t power, igraph_integer_t m, igraph_real_t A);

/**
 * @brief                   This function generates a connected random graph using the Erdos-Renyi model
 * @param n                 The number of vertices in the graph
 * @param type              IGRAPH_ERDOS_RENYI_GNM G(n,m) graph, m edges are selected uniformly randomly in a graph with n vertices\n IGRAPH_ERDOS_RENYI_GNP G(n,p) graph, every possible edge is included in the graph with probability p
 * @param param
 * @return                  Connected random graph using the Erdos-Renyi model
 */
extern igraph_t generateErdosRenyiGraph(igraph_integer_t n, igraph_erdos_renyi_t type, igraph_real_t param);

/**
 * @brief                   This function generates a connected regular random graph
 * @param n                 The number of vertices in the graph
 * @param k                 The degree of each vertex in an undirected graph. For undirected graphs, at least one of k and the number of vertices must be even.
 * @return                  Connected regular random graph
 */
extern igraph_t generateRegularGraph(igraph_integer_t n, igraph_integer_t k);

/**
 * @brief                   This function generate a random graph
 * @param type              Graph distribution:\n 1 geometric\n 2 Barabasi-Albert\n 3 Erdos-Renyi\n 4 regular (clique)
 * @param n                 The number of vertices in the graph
 * @return                  Random graph
 */
extern igraph_t generateRandomGraph(int type, int n);

/**
 * \brief                   This function prints the name of graph distribution.
 * @param type              Graph distribution:\n 1 geometric\n 2 Barabasi-Albert\n 3 Erdos-Renyi\n 4 regular (clique)
 */
extern string printGraphType(int type);
