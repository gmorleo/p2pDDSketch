//
// Created by giuseppe on 16/11/19.
//

#include "graph.h"

int generateGeometricGraph(igraph_t *G_graph, igraph_integer_t n, igraph_real_t radius)
{

    int returnValue = -1;
    igraph_bool_t connected;

    // generate a connected random graph using the geometric model
    returnValue = igraph_grg_game(G_graph, n, radius, 0, nullptr, nullptr);
    if (returnValue){
        printError(GRAPH_GENERATION_ERROR, __FUNCTION__);
        return GRAPH_GENERATION_ERROR;
    }

    igraph_is_connected(G_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(G_graph);
        returnValue = igraph_grg_game(G_graph, n, radius, 0, nullptr, nullptr);
        if (returnValue){
            printError(GRAPH_GENERATION_ERROR, __FUNCTION__);
            return GRAPH_GENERATION_ERROR;
        }

        igraph_is_connected(G_graph, &connected, IGRAPH_WEAK);
    }
    return returnValue;
}

int generateBarabasiAlbertGraph(igraph_t *BA_graph, igraph_integer_t n, igraph_real_t power, igraph_integer_t m, igraph_real_t A)
{

    int returnValue = -1;
    igraph_bool_t connected;

    // generate a connected random graph using the Barabasi-Albert model
    returnValue = igraph_barabasi_game(/* graph=    */ BA_graph,
            /* n=        */ n,
            /* power=    */ power,
            /* m=        */ m,
            /* outseq=   */ 0,
            /* outpref=  */ 0,
            /* A=        */ A,
            /* directed= */ IGRAPH_UNDIRECTED,
            /* algo=     */ IGRAPH_BARABASI_PSUMTREE,
            /* start_from= */ 0);
    if (returnValue){
        printError(GRAPH_GENERATION_ERROR, __FUNCTION__);
        return GRAPH_GENERATION_ERROR;
    }


    igraph_is_connected(BA_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(BA_graph);
        returnValue = igraph_barabasi_game(/* graph=    */ BA_graph,
                /* n=        */ n,
                /* power=    */ power,
                /* m=        */ m,
                /* outseq=   */ 0,
                /* outpref=  */ 0,
                /* A=        */ A,
                /* directed= */ IGRAPH_UNDIRECTED,
                /* algo=     */ IGRAPH_BARABASI_PSUMTREE,
                /* start_from= */ 0);
        if (returnValue){
            printError(GRAPH_GENERATION_ERROR, __FUNCTION__);
            return GRAPH_GENERATION_ERROR;
        }

        igraph_is_connected(BA_graph, &connected, IGRAPH_WEAK);
    }

    return returnValue;
}


int generateErdosRenyiGraph(igraph_t *ER_graph, igraph_integer_t n, igraph_erdos_renyi_t type, igraph_real_t param)
{

    int returnValue = -1;
    igraph_bool_t connected;

    // generate a connected random graph using the Erdos-Renyi model
    returnValue = igraph_erdos_renyi_game(ER_graph, type, n, param, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
    if (returnValue){
        printError(GRAPH_GENERATION_ERROR, __FUNCTION__);
        return GRAPH_GENERATION_ERROR;
    }
    igraph_is_connected(ER_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(ER_graph);
        returnValue = igraph_erdos_renyi_game(ER_graph, type, n, param, IGRAPH_UNDIRECTED, IGRAPH_NO_LOOPS);
        if (returnValue){
            printError(GRAPH_GENERATION_ERROR, __FUNCTION__);
            return GRAPH_GENERATION_ERROR;
        }

        igraph_is_connected(ER_graph, &connected, IGRAPH_WEAK);
    }

    return returnValue;
}

int generateRegularGraph(igraph_t *R_graph, igraph_integer_t n, igraph_integer_t k)
{

    int returnValue = -1;
    igraph_bool_t connected;

    // generate a connected regular random graph
    returnValue = igraph_k_regular_game(R_graph, n, k, IGRAPH_UNDIRECTED, 0);
    if (returnValue){
        printError(GRAPH_GENERATION_ERROR, __FUNCTION__);
        return GRAPH_GENERATION_ERROR;
    }

    igraph_is_connected(R_graph, &connected, IGRAPH_WEAK);
    while(!connected){
        igraph_destroy(R_graph);
        returnValue = igraph_k_regular_game(R_graph, n, k, IGRAPH_UNDIRECTED, 0);
        if (returnValue){
            printError(GRAPH_GENERATION_ERROR, __FUNCTION__);
            return GRAPH_GENERATION_ERROR;
        }

        igraph_is_connected(R_graph, &connected, IGRAPH_WEAK);
    }

    return returnValue;
}

int generateRandomGraph(igraph_t *random_graph, int type, int n)
{
    int returnValue = -1;

    switch (type) {
        case 1:
            returnValue = generateGeometricGraph(random_graph, n, sqrt(100.0/(float)n));
            break;
        case 2:
            returnValue = generateBarabasiAlbertGraph(random_graph, n, 1.0, 5, 1.0);
            break;
        case 3:
            returnValue = generateErdosRenyiGraph(random_graph, n, IGRAPH_ERDOS_RENYI_GNP, 10.0/(float)n);
            // random_graph = generateErdosRenyiGraph(n, IGRAPH_ERDOS_RENYI_GNM, ceil(n^2/3));
            break;
        case 4:
            returnValue = generateRegularGraph(random_graph, n, n-1);
            break;

        default:
            printError(UNKNOWN_GRAPH_TYPE, __FUNCTION__);
            return UNKNOWN_GRAPH_TYPE;
    }

    return returnValue;
}

igraph_t* generateGraph(int peers, int graphType) {

    int returnValue = -1;

    igraph_t* graph = nullptr;
    graph = new (nothrow) igraph_t;
    if (!graph) {
        printError(MEMORY_ERROR,__FUNCTION__);
        return nullptr;
    }

    // turn on attribute handling in igraph
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    // seed igraph PRNG
    igraph_rng_seed(igraph_rng_default(), 42);

    // generate a connected random graph
    returnValue = generateRandomGraph(graph, graphType, peers);
    if (returnValue) {
        return nullptr;
    }

    return graph;
}


string printGraphProperties(igraph_t *graph) {

    igraph_vector_t result;
    igraph_real_t mindeg;
    igraph_real_t maxdeg;
    stringstream string_result;

    igraph_vector_init(&result, 0);
    igraph_degree(graph, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    igraph_vector_minmax(&result, &mindeg, &maxdeg);

    string_result << "Graph properties:\nminimum degree = "  <<  mindeg << ";\nmaximum degree = " << maxdeg << ";\n\n";

    igraph_vector_destroy(&result);

    return string_result.str();
}

string printGraphType(int type)
{
    string graph_name;
    switch (type) {
        case 1:
            graph_name = "Geometric random graph\n";
        case 2:
            graph_name = "Barabasi-Albert random graph\n";
        case 3:
            graph_name = "Erdos-Renyi random graph\n";
        case 4:
            graph_name = "Regular random graph\n";
    }

    return graph_name;
}
