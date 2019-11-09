/********************************************************************
 DDSketch

 An algorithm for tracking quantiles in data streams

 Charles Masson, Jee E. Rim, and Homin K. Lee. 2019. DDSketch: a fast and fully-mergeable quantile sketch with relative-error guarantees. Proc. VLDB Endow. 12, 12 (August 2019), 2195-2205. DOI: https://doi.org/10.14778/3352063.3352135

 This implementation by
 by Giuseppe Morleo
 University of Salento, Italy

 *********************************************************************/
/// \file


#include <iostream>
#include <igraph/igraph.h>
#include <cstring>
#include <random>
#include <map>
#include <chrono>
#include <fstream>
#include "ddsketch.h"

using namespace std;
using namespace std::chrono;

const int DEFAULT_OFFSET = 1073741824; //2^31/2
const int DEFAULT_BIN_LIMIT = 500;
const float DEFAULT_ALPHA = 0.01;

void startTheClock();
double stopTheClock();
void usage(char* cmd);
igraph_t generateGeometricGraph(igraph_integer_t n, igraph_real_t radius);
igraph_t generateBarabasiAlbertGraph(igraph_integer_t n, igraph_real_t power, igraph_integer_t m, igraph_real_t A);
igraph_t generateErdosRenyiGraph(igraph_integer_t n, igraph_erdos_renyi_t type, igraph_real_t param);
igraph_t generateRegularGraph(igraph_integer_t n, igraph_integer_t k);
igraph_t generateRandomGraph(int type, int n);
void printGraphType(int type);

long getDatasetSize(const string &name_file);
int loadDataset(const string &name_file, double *dataset);

struct Params {
    /// Number of points in the dataset
    long        ni;
    /// Name of dataset
    string      name_file;
    /// Number of possible distinct items
    uint32_t    domainSize;
    /// Graph distribution: 1 geometric 2 Barabasi-Albert 3 Erdos-Renyi 4 regular (clique)
    int         graphType;
    /// Number of peer in the net
    int         peers;
    /// Fan-out of peers
    int         fanOut;
    /// Threshold to check the peer's convergence
    double      convThreshold;
    /// Number of consecutive rounds in which a peer must locally converge
    int         convLimit;
    /// Fixed number of round to execute
    int         roundsToExecute;
    /// Output file to redirect stdout
    string      outputFilename;
};

high_resolution_clock::time_point t1, t2;

int main(int argc, char **argv) {
    /*** Default Parameters***/

    long        ni;                         // points number
    string      name_file = "../normal_mean_2_stdev_3.csv";
    uint32_t    domainSize = 1048575;       // number of possible distinct items
    int         graphType = 2;              // graph distribution: 1 geometric 2 Barabasi-Albert 3 Erdos-Renyi 4 regular (clique)
    int         peers = 1000;                 // number of peers
    int         fanOut = 5;                 // fan-out of peers
    double      convThreshold = 0.0001;     // local convergence tolerance
    int         convLimit = 3;              // number of consecutive rounds in which a peer must locally converge
    int         roundsToExecute = -1;       // Se metto un numero esegue quello

    float       alpha = DEFAULT_ALPHA;
    int         offset = DEFAULT_OFFSET;
    int         bin_limit = DEFAULT_BIN_LIMIT;

    int         p_star = -1;
    double      delta = 0.04;
    double          elapsed;
    int             iterations;


    bool            autoseed = false;
    bool        outputOnFile = false;
    string      outputFilename;
    long        *peerLastItem;              // index of a peer last item
    igraph_t    graph;
    Params      params;

    /*** parse command-line parameters ***/
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-di") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing domain size parameter." << endl;
                return -1;
            }
            domainSize = stol(argv[i]);
        } else if (strcmp(argv[i], "-delta") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing delta parameter." << endl;
                return -1;
            }
            delta = stod(argv[i]);
        } else if (strcmp(argv[i], "-p") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of peers parameter." << endl;
                return -1;
            }
            peers = stoi(argv[i]);
        } else if (strcmp(argv[i], "-ps") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing maximum number of peers parameter." << endl;
                return -1;
            }
            p_star = stoi(argv[i]);
        } else if (strcmp(argv[i], "-f") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing fan-out parameter." << endl;
                return -1;
            }
            fanOut = stoi(argv[i]);;
        } else if (strcmp(argv[i], "-d") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing graph type parameter" << endl;
                return -1;
            }
            graphType = stoi(argv[i]);
        } else if (strcmp(argv[i], "-ct") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing convergence tolerance parameter." << endl;
                return -1;
            }
            convThreshold = stod(argv[i]);
        } else if (strcmp(argv[i], "-cl") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing # of consecutive rounds in which convergence is satisfied parameter." << endl;
                return -1;
            }
            convLimit = stol(argv[i]);
        } else if (strcmp(argv[i], "-of") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing filename for simulation output." << endl;
                return -1;
            }
            outputFilename = string(argv[i]);
        } else if (strcmp(argv[i], "-r") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of rounds to execute.\n";
                return -1;
            }
            roundsToExecute = stoi(argv[i]);
        } else if (strcmp(argv[i], "-al") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of accuracy for DDSketch.\n";
                return -1;
            }
            alpha = stod(argv[i]);
        } else if (strcmp(argv[i], "-ofs") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of offset for DDSketch.\n";
                return -1;
            }
            offset = stoi(argv[i]);
        } else if (strcmp(argv[i], "-bl") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of bins limit for DDSketch.\n";
                return -1;
            }
            bin_limit = stoi(argv[i]);
        } else if (strcmp(argv[i], "-dt") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing file name of dataset.\n";
                return -1;
            }
            name_file = argv[i];
        }else if (strcmp(argv[i], "-as") == 0) {
            autoseed = true;
        } else {
            //usage(argv[0]);
            return -1;
        }
    }

    /*** read dataset dimensions ***/
    ni = getDatasetSize(name_file);
    double *dataset;
    dataset = new double[ni];
    loadDataset(name_file, dataset);

    /*** Compute last item for each peer***/
    peerLastItem = new long[peers]();
    std::random_device rd;                                   // obtain a random number from hardware
    std::mt19937 eng(rd());                                  // seed the generator
    std::uniform_real_distribution<> distr(-1, 1);    // define the range

    for(int i = 0; i < peers - 1; i++){
        float rnd = distr(eng);
        //cerr << "rnd: " << rnd << "\n";
        long last_item = rnd * ((float)ni/(float)peers) * 0.1 + (float) (i+1) * ((float)ni/(float)peers) - 1;
        peerLastItem[i] = last_item;
    }

    peerLastItem[peers - 1] = ni-1;

    /*** check the partitioning correctness ***/
    int lastNotNull=0;
    long sum = peerLastItem[0] + 1;
    //cerr << "peer 0:" << sum << "\n";
    for(int i = 1; i < peers; i++) {
        if (!peerLastItem[i] && peerLastItem[i] != peerLastItem[i-1]) {
            lastNotNull = i-1;
            //cerr << "peer " << i << ":" << 0 << "\n";
        } else {
            if (peerLastItem[i-1] != 0) {
                sum += peerLastItem[i] - peerLastItem[i-1];
                //cerr << "peer " << i << ":" << peerLastItem[i] - peerLastItem[i-1] << "\n";
            }
            else if (peerLastItem[i]){
                sum += peerLastItem[i] - peerLastItem[lastNotNull];
                //cerr << "peer " << lastNotNull+1 << ":" << peerLastItem[i] - peerLastItem[lastNotNull] << "\n";
            }
        }
    }
    //cout<< "sum: " << sum << endl;
    if(sum != ni) {
        cerr << "ERROR: ni = "<< ni << "!= sum = " << sum << endl;
        exit(EXIT_FAILURE);
    }
    cout.flush();

    /*** assign parameters read from command line ***/
    params.ni = ni;
    params.name_file = name_file;
    params.domainSize = domainSize;
    params.graphType = graphType;
    params.peers = peers;
    params.fanOut = fanOut;
    params.convThreshold = convThreshold;
    params.convLimit = convLimit;
    params.roundsToExecute = roundsToExecute;
    params.outputFilename = outputFilename;

    outputOnFile = !params.outputFilename.empty();
    if (!outputOnFile) {
        printf("\n\nPARAMETERS:\n");
        //printf("distinct items in the stream = %d\n", params.domainSize);
        cout << "dataset = " << params.name_file << "\n";
        cout << "n° points = " << params.ni << "\n";
        cout << "alpha = " << alpha << "\n";
        cout << "offset = " << offset << "\n";
        cout << "bin_limit = " << bin_limit << "\n";
        cout << "peers = " << params.peers << "\n";
        cout << "fan-out = " << params.fanOut << "\n";
        cout << "graph type = ";
        printGraphType(params.graphType);
        cout << "local convergence tolerance = "<< params.convThreshold << "\n";
        cout << "number of consecutive rounds in which a peer must locally converge = "<< params.convLimit << "\n";
        cout << "\n\n";
    }

    /*** Graph generation ***/
    // turn on attribute handling in igraph
    igraph_i_set_attribute_table(&igraph_cattribute_table);

    // seed igraph PRNG
    igraph_rng_seed(igraph_rng_default(), 42);

    startTheClock();
    // generate a connected random graph
    graph = generateRandomGraph(params.graphType, params.peers);

    double graphgentime = stopTheClock();
    if (!outputOnFile) {
        cout <<"Time (seconds) required to generate the random graph: " << graphgentime << "\n";
    }

    // determine minimum and maximum vertex degree for the graph
    igraph_vector_t result;
    igraph_real_t mindeg;
    igraph_real_t maxdeg;

    igraph_vector_init(&result, 0);
    igraph_degree(&graph, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);
    igraph_vector_minmax(&result, &mindeg, &maxdeg);
    if (!outputOnFile) {
        cout << "Minimum degree is "  <<  mindeg << ", Maximum degree is " << maxdeg << "\n";
    }

    /*** Declaration of peer sketches ***/
    auto *dds = new DDS_type*[params.peers];

    startTheClock();
    /*** Distributed computation simulation ***/
    long start = 0;
    for(int peerID = 0; peerID < params.peers; peerID++){

        dds[peerID] = DDS_Init(offset, bin_limit, alpha);

        for ( long i = start; i <= peerLastItem[peerID]; i++ ) {
            DDS_Add(dds[peerID], dataset[i]);
        }

        //cout << "PeerID = " << peerID << " sketch size = " << DDS_Size(dds[peerID]) << " alpha = " << dds[peerID]->alpha << endl;
        start = peerLastItem[peerID] + 1;
    }

    double distributed_time = stopTheClock();
    if (!outputOnFile) {
        cout <<"Time (seconds) required to add all elements for all the peers: " << distributed_time << "\n";
    }

    // Weight for sum only one peer is setted to 1
    auto *dimestimate = new double[params.peers]();
    dimestimate[0] = 1;

    // Weight at rount t-1
    auto *prevestimate = new double[params.peers]();

    // Number of peers that have not reached convergence
    int Numberofconverged = params.peers;

    // Converged peers
    auto *converged = new bool[params.peers]();
    for(int i = 0; i < params.peers; i++)
        converged[i] = false;

    // Local convergence tolerance
    auto *convRounds = new int[params.peers]();

    int rounds = 0;

    if (!outputOnFile) {
        cout <<"\nStarting distributed agglomeration merge..." << endl;
    }

    startTheClock();
    /*** Merge information about agglomeration ***/
    while( (params.roundsToExecute < 0 && Numberofconverged) || params.roundsToExecute > 0) {

        memcpy(prevestimate, dimestimate, params.peers * sizeof(double));

        for (int peerID = 0; peerID < params.peers; peerID++) {

            // determine peer neighbors
            igraph_vector_t neighbors;
            igraph_vector_init(&neighbors, 0);
            igraph_neighbors(&graph, &neighbors, peerID, IGRAPH_ALL);
            long neighborsSize = igraph_vector_size(&neighbors);
            if (fanOut < neighborsSize && fanOut != -1) {
                // randomly sample f adjacent vertices
                igraph_vector_shuffle(&neighbors);
                igraph_vector_remove_section(&neighbors, params.fanOut, neighborsSize);
            }

            neighborsSize = igraph_vector_size(&neighbors);
            for (long i = 0; i < neighborsSize; i++) {
                int neighborID = (int) VECTOR(neighbors)[i];
                igraph_integer_t edgeID;
                igraph_get_eid(&graph, &edgeID, peerID, neighborID, IGRAPH_UNDIRECTED, 1);

                //DDS_mergeGossip(dds[peerID], dds[neighborID]);
                DDS_merge(dds[peerID], dds[neighborID]);
                DDS_replaceBinMap(dds[peerID], dds[neighborID]);

                double mean = (dimestimate[peerID] + dimestimate[neighborID]) / 2;
                dimestimate[peerID] = mean;
                dimestimate[neighborID] = mean;
            }

            igraph_vector_destroy(&neighbors);
        }

        // check local convergence, if roundsToExecute is less than 0, the algorithm will be running until convergence
        if (params.roundsToExecute < 0) {

            for(int peerID = 0; peerID < params.peers; peerID++){

                if(converged[peerID])
                    continue;

                // Check local convergence
                bool dimestimateconv;
                if(prevestimate[peerID])
                    dimestimateconv = fabs((prevestimate[peerID] - dimestimate[peerID]) / prevestimate[peerID]) < params.convThreshold;
                else
                    dimestimateconv = false;

                // Increase rounds of convergence
                if(dimestimateconv)
                    convRounds[peerID]++;
                else
                    convRounds[peerID] = 0;

                //printf ("PeerID %d, round %d, convRound %d\n", peerID, rounds, convRounds[peerID]);

                // If a peer reaches convergence, decrease by one the numberofconverged
                converged[peerID] = (convRounds[peerID] >= params.convLimit);
                if(converged[peerID]){
                    //printf("peer %d rounds before convergence: %d\n", peerID, rounds + 1);
                    Numberofconverged --;
                }
            }
        }
        rounds++;
        cerr << "\r Active peers: " << Numberofconverged << " - Rounds: " << rounds << "          " << endl;
        params.roundsToExecute--;
    }
    double comunication_time = stopTheClock();
    if (!outputOnFile) {
        cout <<"Time (seconds) required to reach convergence: " << comunication_time << "\n";
    }

    cout << "We are: " << 1/dimestimate[0] << endl;

    /*** Distributed computation simulation ***/
    // Finalize the sum by dividing each value by dimestimate[peerID]
    for(int peerID = 0; peerID < params.peers; peerID++){
        DDS_finalizeGossip(dds[peerID], dimestimate[peerID]);
    }

    // determine the 0.7 quantile
    float q = 0.7;
    int idx = floor(1+q*(params.ni-1));

    // determine the correct answer
    // i.e., the number stored at the index (idx-1) in the sorted permutation of the vector
    // note that we are not sorting the vector, we are using the quickselect() algorithm
    // which in C++ is available as std::nth_element
    nth_element(dataset, dataset + (idx-1), dataset +  params.ni);

    double quantile = DDS_GetQuantile(dds[0], q);
    double error = abs((quantile-dataset[idx-1])/dataset[idx-1]);

    cout << "Result for q = " << q << endl;
    cout << "Real: " << dataset[idx-1] << " Estimation: " << quantile << " Error: " << error << endl;


    return 0;
}

/**
 * @brief                   This function sets the start time
 */
void startTheClock(){
    t1 = high_resolution_clock::now();
}

/**
 * @brief                   This function returns the time between the start time and the end time
 * @return                  Total time between two times
 */
double stopTheClock() {
    t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    return time_span.count();
}

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

/**
 * \brief                   This function prints the name of graph distribution.
 * @param type              Graph distribution:\n 1 geometric\n 2 Barabasi-Albert\n 3 Erdos-Renyi\n 4 regular (clique)
 */
void printGraphType(int type)
{

    switch (type) {
        case 1:
            cout << "Geometric random graph\n";
            break;
        case 2:
            cout << "Barabasi-Albert random graph\n";
            break;
        case 3:
            cout << "Erdos-Renyi random graph\n";
            break;
        case 4:
            cout << "Regular random graph\n";
            break;

        default:
            break;
    }

}

/**
 * \brief                   This function computes the dimension of the dataset
 * @param name_file         Name of dataset
 * @return                  Return the number of element in the dataset(row)
 */
long getDatasetSize(const string &name_file) {

    ifstream inputFile(name_file);
    string line;

    long rows = 0;

    while (getline(inputFile, line))
        rows++;

    return rows;
}

/**
 * \brief                   This function loads the dataset into an array
 * @param name_file         Name of dataset
 * @param dataset           Array
 * @return                  An array containing the whole dataset
 */
int loadDataset(const string &name_file, double *dataset) {

    ifstream inputFile(name_file);
    string line;

    int row = 0;

    while (getline(inputFile, line)) {
        dataset[row] = stod(line);
        row++;
    }

    return 0;
}