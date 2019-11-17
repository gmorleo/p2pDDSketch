//
// Created by giuseppe on 16/11/19.
//

#include "graph.h"

/**
 * @brief                   This function generates a connetected random graph using the geometric model
 * @param n                 The number of vertices in the graph
 * @param radius            Graph radius
 * @return                  Connected random graph using the geometric model
 */
igraph_t generateGeometricGraph(igraph_integer_t n, igraph_real_t radius)
{
    igraph_t G_graph;
    igraph_bool_t connected;

    // generate a connected random graph using the geometric model
    igraph_grg_game(&G_graph, n, radius, 0, nullptr, nullptr);

    igraph_is_connected(&G_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(&G_graph);
        igraph_grg_game(&G_graph, n, radius, 0, nullptr, nullptr);

        igraph_is_connected(&G_graph, &connected, IGRAPH_WEAK);
    }

    return G_graph;
}

/**
 * @brief                   This function generates a connected random graph using the Barabasi-Albert model
 * @param n                 The number of vertices in the graph
 * @param power             Power of the preferential attachment. The probability that a vertex is cited is proportional to d^power+A, where d is its degree, power and A are given by arguments. In the classic preferential attachment model power=1
 * @param m                 m = number of outgoing edges generated for each vertex
 * @param A                 A = The probability that a vertex is cited is proportional to d^power+A, where d is its degree, power and A are given by argument
 * @return                  Connected random graph using the Barabasi-Albert model
 */
igraph_t generateBarabasiAlbertGraph(igraph_integer_t n, igraph_real_t power, igraph_integer_t m, igraph_real_t A)
{

    igraph_t BA_graph;
    igraph_bool_t connected;

    // generate a connected random graph using the Barabasi-Albert model
    igraph_barabasi_game(/* graph=    */ &BA_graph,
            /* n=        */ n,
            /* power=    */ power,
            /* m=        */ m,
            /* outseq=   */ 0,
            /* outpref=  */ 0,
            /* A=        */ A,
            /* directed= */ IGRAPH_UNDIRECTED,
            /* algo=     */ IGRAPH_BARABASI_PSUMTREE,
            /* start_from= */ 0);


    igraph_is_connected(&BA_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(&BA_graph);
        igraph_barabasi_game(/* graph=    */ &BA_graph,
                /* n=        */ n,
                /* power=    */ power,
                /* m=        */ m,
                /* outseq=   */ 0,
                /* outpref=  */ 0,
                /* A=        */ A,
                /* directed= */ IGRAPH_UNDIRECTED,
                /* algo=     */ IGRAPH_BARABASI_PSUMTREE,
                /* start_from= */ 0);

        igraph_is_connected(&BA_graph, &connected, IGRAPH_WEAK);
    }

    return BA_graph;
}

/**
 * @brief                   This function generates a connected random graph using the Erdos-Renyi model
 * @param n                 The number of vertices in the graph
 * @param type              IGRAPH_ERDOS_RENYI_GNM G(n,m) graph, m edges are selected uniformly randomly in a graph with n vertices\n IGRAPH_ERDOS_RENYI_GNP G(n,p) graph, every possible edge is included in the graph with probability p
 * @param param
 * @return                  Connected random graph using the Erdos-Renyi model
 */
igraph_t generateErdosRenyiGraph(igraph_integer_t n, igraph_erdos_renyi_t type, igraph_real_t param)
{

    igraph_t ER_graph;
    igraph_bool_t connected;

    // generate a connected random graph using the Erdos-Renyi model
    igraph_erdos_renyi_game(&ER_graph, type, n, param, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

    igraph_is_connected(&ER_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(&ER_graph);
        igraph_erdos_renyi_game(&ER_graph, type, n, param, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);

        igraph_is_connected(&ER_graph, &connected, IGRAPH_WEAK);
    }

    return ER_graph;
}

/**
 * @brief                   This function generates a connected regular random graph
 * @param n                 The number of vertices in the graph
 * @param k                 The degree of each vertex in an undirected graph. For undirected graphs, at least one of k and the number of vertices must be even.
 * @return                  Connected regular random graph
 */
igraph_t generateRegularGraph(igraph_integer_t n, igraph_integer_t k)
{

    igraph_t R_graph;
    igraph_bool_t connected;

    // generate a connected regular random graph
    igraph_k_regular_game(&R_graph, n, k, IGRAPH_UNDIRECTED, 0);

    igraph_is_connected(&R_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(&R_graph);
        igraph_k_regular_game(&R_graph, n, k, IGRAPH_UNDIRECTED, 0);

        igraph_is_connected(&R_graph, &connected, IGRAPH_WEAK);
    }

    return R_graph;
}

/**
 * @brief                   This function generate a random graph
 * @param type              Graph distribution:\n 1 geometric\n 2 Barabasi-Albert\n 3 Erdos-Renyi\n 4 regular (clique)
 * @param n                 The number of vertices in the graph
 * @return                  Random graph
 */
igraph_t generateRandomGraph(int type, int n)
{
    igraph_t random_graph;

    switch (type) {
        case 1:
            random_graph = generateGeometricGraph(n, sqrt(100.0/(float)n));
            break;
        case 2:
            random_graph = generateBarabasiAlbertGraph(n, 1.0, 5, 1.0);
            break;
        case 3:
            random_graph = generateErdosRenyiGraph(n, IGRAPH_ERDOS_RENYI_GNP, 10.0/(float)n);
            // random_graph = generateErdosRenyiGraph(n, IGRAPH_ERDOS_RENYI_GNM, ceil(n^2/3));
            break;
        case 4:
            random_graph = generateRegularGraph(n, n-1);
            break;

        default:
            break;
    }

    return random_graph;

}
/*

*/
/**
 * \brief                   This function prints the name of graph distribution.
 * @param type              Graph distribution:\n 1 geometric\n 2 Barabasi-Albert\n 3 Erdos-Renyi\n 4 regular (clique)
 */

string printGraphType(int type)
{

    switch (type) {
        case 1:
            return "Geometric random graph\n";
        case 2:
            return "Barabasi-Albert random graph\n";
        case 3:
            return "Erdos-Renyi random graph\n";
        case 4:
            return "Regular random graph\n";

        default:
            return "\n";
    }

}
