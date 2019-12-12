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
#include <chrono>
#include <fstream>
#include <iomanip>
#include <codecvt>
#include "ddsketch.h"
#include "graph.h"

#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define GREEN       "\x1B[32m"
#define RESET       "\033[0m"

using namespace std;
using namespace std::chrono;

const string        DEFAULT_FILENAME = "../normal_mean_2_stdev_3.csv";
const string        DEFAULT_OUTPUT_FILE = "out.txt";
const uint32_t      DEFAULT_DOMAIN_SIZE = 1048575;
const int           DEFAULT_GRAPH_TYPE = 2;
const int           DEFAULT_PEERS = 1000;
const int           DEFAULT_FAN_OUT = 5;
const double        DEFAULT_CONVERGENCE_THRESHOLD = 0.0001;
const int           DEFAULT_CONVERGENCE_LIMIT = 3;
const int           DEFAULT_ROUND_TO_EXECUTE = -1;
const bool          DEFAULT_AUTO_SEED = false;
const int           DEFAULT_OFFSET = 1073741824; //2^31/2
const int           DEFAULT_BIN_LIMIT = 500;
const float         DEFAULT_ALPHA = 0.008;

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
    ///
    bool        outputOnFile;
    /// Output file to redirect stdout
    string      outputFilename;
    ///
    bool        autoseed;

    /// used to allow the skecth storing both positive, 0 and negative values
    int         offset;
    /// maximum number of bins
    int         bin_limit;
    /// this parameter defines alpha-accuracy of a q-quantile
    double      alpha;
};

/**
 * @brief                   This function sets the start time
 */
void startTheClock();

/**
 * @brief                   This function returns the time between the start time and the end time
 * @return                  Total time between two times
 */
double stopTheClock();

/**
 * @brief                   Parmas costructor
 * @return                  An allocated Params data structure with default parameters
 */
Params* init();

/**
 * @brief                   This function parses command-line arguments
 * @param argc              Count of command line arguments
 * @param argv              Command line arguments
 * @param params            Params data structure
 * @return                  0 in case of success, -13 in case of usage error
 */
int parse(int argc, char** argv, Params* params);

/**
 * @brief                   This function prints the optional parameters of the application
 * @param cmd               The name of the application.
 */
void usage(char* cmd);

/**
 * \brief                   This function computes the dimension of the dataset
 * @param name_file         Name of dataset
 * @param row               Where the number of row is stored
 * @return                  0 in case of success, -3 in case of file error
 */
int getDatasetSize(const string &name_file, int &rows);

/**
 * \brief                   This function loads the dataset into an array
 * @param name_file         Name of dataset
 * @param dataset           Array
 * @return                  An array containing the whole dataset
 */
int loadDataset(const string &name_file, double *dataset);

/**
 * @brief
 * @param params
 * @param peerLastItem
 * @return
 */
int computeLastItem(Params* params, long* peerLastItem);

/**
 * @brief
 * @param params
 * @return
 */
int printParameters(Params* params);

string getParameters(Params* params);

/**
 * @brief
 * @param params
 * @param dds
 * @return
 */
int distributedInitializeSketch(Params* params, DDS_type** dds);

/**
 * @brief
 * @param params
 * @param dds
 * @param dataset
 * @param peerLastItem
 * @return
 */
int distributedAdd(Params* params, DDS_type** dds, double* dataset, long* peerLastItem);

/**
 * @brief
 * @param params
 * @param dds
 * @param graph
 * @return
 */
int distributedCommunication(Params* params, DDS_type** dds, igraph_t graph);

/**
 * @brief
 * @param params
 * @param dds
 * @param dimestimate
 * @return
 */
int distributedFinalizeMerge(Params* params, DDS_type** dds, double* dimestimate);

/**
 * @brief               This function computes the quantile
 * @param dds           Parameters of the sketch
 * @param stream        Vector that contains all the real values inserted
 * @param n_element     Number of element
 * @return              0: success; \n-2: error;
 */
int testQuantile(DDS_type *dds, double* stream, long n_element, stringstream &output);

high_resolution_clock::time_point t1, t2;

int main(int argc, char **argv) {

    // Declaration

    int returnValue = -1;

    ofstream out;
    stringstream output;

    Params *params = nullptr;
    DDS_type** dds = nullptr;

    double* dataset = nullptr;
    long* peerLastItem = nullptr;

    igraph_t graph;

    /*** Initialize the alghoritm based on default parameters ***/
    params = init();
    if (!params) {
        printError(PARAM_ERROR, __FUNCTION__);
        returnValue = PARAM_ERROR;
        goto ON_EXIT;
    }

    /*** Parse user-supplied parameters ***/
    if ( argc > 1) {
        returnValue = parse(argc, argv, params);
        if (returnValue) {
            goto ON_EXIT;
        }
    }

    params->outputFilename = DEFAULT_OUTPUT_FILE;
    params->outputOnFile = true;

    /*** Open file for output ***/
    if (params->outputOnFile) {
        out.open(params->outputFilename);
        if (!out.is_open()) {
            printError(FILE_ERROR, __FUNCTION__);
            returnValue = FILE_ERROR;
            goto ON_EXIT;
        }
    }

    if(!params->outputOnFile) {
        cout << BOLDRED << "\nPARAMETERS:\n" << RESET;
        cout << getParameters(params);
    } else {
        out << getParameters(params);
    }

    /*** Initilize dataset array ***/
    dataset = new (nothrow) double[params->ni];
    if (!dataset) {
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }
    returnValue = loadDataset(params->name_file, dataset);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Compute last item for each peer***/
    peerLastItem = new (nothrow) long[params->peers]();
    if (!peerLastItem) {
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;;
    }

    returnValue = computeLastItem(params, peerLastItem);
    if (returnValue) {
        goto ON_EXIT;
    }

    /*** Generate random Graph ***/
    returnValue = generateGraph(graph, params->peers, params->graphType);
    if (returnValue) {
        goto ON_EXIT;
    }

    if(!params->outputOnFile) {
        cout << printGraphProperties(graph);
    } else {
        out << printGraphProperties(graph);
    }

    /*** Declaration of peer sketches ***/
    dds = new (nothrow) DDS_type* [params->peers];
    if (!dds) {
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }


    /*** Distributed computation ***/
    returnValue = distributedInitializeSketch(params, dds);
    if (returnValue) {
        cout << "Error with the add" << endl;
        goto ON_EXIT;
    }

    /*** Distributed computation ***/
    returnValue = distributedAdd(params, dds, dataset, peerLastItem);
    if (returnValue) {
        cout << "Error with the add" << endl;
        goto ON_EXIT;
    }

    /*** Distributed communication ***/
    returnValue = distributedCommunication(params, dds, graph);
    if ( returnValue < 0 ) {
        cout << "Error during the distributed communication" << endl;
        goto ON_EXIT;
    }

    /*** Quantile computing ***/
    output.str("");
    returnValue = testQuantile(dds[0], dataset, params->ni, output);
    if ( returnValue < 0 ) {
        cout << "Error during the computation of quantile" << endl;
        goto ON_EXIT;
    }

    if(!params->outputOnFile) {
        cout << output.str();
    } else {
        out << output.str();
    }

    ON_EXIT:

    if(dataset != nullptr) {
        delete[] dataset, dataset = nullptr;
    }

    if(peerLastItem != nullptr) {
        delete[] peerLastItem, peerLastItem = nullptr;
    }

    if(dds != nullptr) {
        for (int i = 0; i < params->peers; ++i) {
            DDS_Destroy(dds[i]), dds[i] = nullptr;
        }
        delete dds, dds = nullptr;
    }

    if(params != nullptr) {
        delete  params, params = nullptr;
    }

    if (out.is_open()) {
        out.close();
    }

    igraph_destroy(&graph);

    return returnValue;
}

int distributedInitializeSketch(Params* params, DDS_type** dds) {

    cout << BOLDRED << "\nStart distributed inizialization.." << RESET << endl;

    int returnValue = -1;

    startTheClock();

    for(int peerID = 0; peerID < params->peers; peerID++) {

        // Initialization peers sketches
        dds[peerID] = DDS_Init(params->offset, params->bin_limit, params->alpha);
        if (!dds[peerID]) {
            return MEMORY_ERROR;
        }

    }

    double distributed_time = stopTheClock();
    cout <<"Time (seconds) required to initialize all sketch for all the peers: " << distributed_time << "\n";

    return  SUCCESS;
}



int distributedAdd(Params* params, DDS_type** dds, double* dataset, long* peerLastItem) {

    cout << BOLDRED << "\nStart distributed computation.." << RESET << endl;

    int returnValue = -1;

    startTheClock();
    long start = 0;
    for(int peerID = 0; peerID < params->peers; peerID++){

        // Adding elements to the sketch
        for ( long i = start; i <= peerLastItem[peerID]; i++ ) {
            returnValue = DDS_AddCollapse(dds[peerID], dataset[i]);
            if (returnValue) {
                return returnValue;
            }
        }

        start = peerLastItem[peerID] + 1;
    }

    double distributed_time = stopTheClock();
    cout <<"Time (seconds) required to add all elements for all the peers: " << distributed_time << "\n";

    return returnValue;
}

int distributedFinalizeMerge(Params* params, DDS_type** dds, double* dimestimate) {

    int returnValue = -1;

    // Finalize the sum by dividing each value by dimestimate[peerID]
    for(int peerID = 0; peerID < params->peers; peerID++){
        returnValue = DDS_finalizeMerge(dds[peerID], dimestimate[peerID]);
        if (returnValue) {
            return returnValue;
        }
    }

    return returnValue;
}

int distributedCommunication(Params* params, DDS_type** dds, igraph_t graph) {

    cout << BOLDRED << "\nStart distributed communication.." << RESET << endl;

    int returnValue = -1;

    double* dimestimate = nullptr;
    double* prevestimate = nullptr;
    bool* converged = nullptr;
    int* convRounds = nullptr;

    double comunication_time;
    int rounds = 0;
    // Number of peers that have not reached convergence
    int Numberofconverged = params->peers;

    // Weight for sum only one peer is setted to 1
    dimestimate = new (nothrow) double[params->peers]();
    if(!dimestimate) {
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }
    dimestimate[0] = 1;

    // Weight at rount t-1
    prevestimate = new (nothrow) double[params->peers]();
    if(!prevestimate) {
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }

    // Converged peers
    converged = new (nothrow) bool[params->peers]();
    if(!converged) {
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }

    for(int i = 0; i < params->peers; i++) {
        converged[i] = false;
    }

    // Local convergence tolerance
    convRounds = new (nothrow) int[params->peers]();
    if(!convRounds) {
        printError(MEMORY_ERROR, __FUNCTION__);
        returnValue = MEMORY_ERROR;
        goto ON_EXIT;
    }

    if (!params->outputOnFile) {
        cout <<"Starting distributed agglomeration merge..." << endl;
    }

    startTheClock();
    /*** Merge information about agglomeration ***/
    while( (params->roundsToExecute < 0 && Numberofconverged) || params->roundsToExecute > 0) {

        memcpy(prevestimate, dimestimate, params->peers * sizeof(double));

        for (int peerID = 0; peerID < params->peers; peerID++) {

            // determine peer neighbors
            igraph_vector_t neighbors;
            igraph_vector_init(&neighbors, 0);
            igraph_neighbors(&graph, &neighbors, peerID, IGRAPH_ALL);
            long neighborsSize = igraph_vector_size(&neighbors);
            if (params->fanOut < neighborsSize && params->fanOut != -1) {
                // randomly sample f adjacent vertices
                igraph_vector_shuffle(&neighbors);
                igraph_vector_remove_section(&neighbors, params->fanOut, neighborsSize);
            }

            neighborsSize = igraph_vector_size(&neighbors);
            for (long i = 0; i < neighborsSize; i++) {
                int neighborID = (int) VECTOR(neighbors)[i];
                igraph_integer_t edgeID;
                igraph_get_eid(&graph, &edgeID, peerID, neighborID, IGRAPH_UNDIRECTED, 1);

                returnValue = DDS_MergeCollapse(dds[peerID], dds[neighborID]);
                if (returnValue) {
                    goto ON_EXIT;
                }

                returnValue = DDS_replaceBinMap(dds[peerID], dds[neighborID]);
                if (returnValue) {
                    goto ON_EXIT;
                }

                double mean = (dimestimate[peerID] + dimestimate[neighborID]) / 2;
                dimestimate[peerID] = mean;
                dimestimate[neighborID] = mean;
            }

            igraph_vector_destroy(&neighbors);
        }

        // check local convergence, if roundsToExecute is less than 0, the algorithm will be running until convergence
        if (params->roundsToExecute < 0) {

            for(int peerID = 0; peerID < params->peers; peerID++){

                if(converged[peerID])
                    continue;

                // Check local convergence
                bool dimestimateconv;
                if(prevestimate[peerID])
                    dimestimateconv = fabs((prevestimate[peerID] - dimestimate[peerID]) / prevestimate[peerID]) < params->convThreshold;
                else
                    dimestimateconv = false;

                // Increase rounds of convergence
                if(dimestimateconv)
                    convRounds[peerID]++;
                else
                    convRounds[peerID] = 0;

                //printf ("PeerID %d, round %d, convRound %d\n", peerID, rounds, convRounds[peerID]);

                // If a peer reaches convergence, decrease by one the numberofconverged
                converged[peerID] = (convRounds[peerID] >= params->convLimit);
                if(converged[peerID]){
                    //printf("peer %d rounds before convergence: %d\n", peerID, rounds + 1);
                    Numberofconverged --;
                }
            }
        }
        rounds++;
        cout << GREEN <<" Active peers: " << Numberofconverged << " - Rounds: " << rounds << RESET << endl;
        params->roundsToExecute--;
    }
    comunication_time = stopTheClock();
    cout <<"\nTime (seconds) required to reach convergence: " << comunication_time << "\n";

    returnValue = distributedFinalizeMerge(params, dds, dimestimate);
    if (returnValue) {
        goto ON_EXIT;
    }

    cout << "Estimate of peers number: " << 1/dimestimate[0] << endl;

    ON_EXIT:

    if (dimestimate != nullptr) {
        delete[] dimestimate, dimestimate = nullptr;
    }
    if ( prevestimate != nullptr ) {
        delete[] prevestimate, prevestimate = nullptr;
    }
    if ( converged != nullptr ) {
        delete[] converged, converged = nullptr;
    }
    if ( convRounds != nullptr ) {
        delete[] convRounds, convRounds = nullptr;
    }

    return returnValue;
}

int computeLastItem(Params* params, long* peerLastItem) {

    std::random_device rd;                                   // obtain a random number from hardware
    std::mt19937 eng(rd());                                  // seed the generator
    std::uniform_real_distribution<> distr(-1, 1);    // define the range

    for(int i = 0; i < params->peers - 1; i++){
        float rnd = distr(eng);
        //cerr << "rnd: " << rnd << "\n";
        long last_item = rnd * ((float)params->ni/(float)params->peers) * 0.1 + (float) (i+1) * ((float)params->ni/(float)params->peers) - 1;
        peerLastItem[i] = last_item;
    }

    peerLastItem[params->peers - 1] = params->ni-1;

    /*** check the partitioning correctness ***/
    int lastNotNull=0;
    long sum = peerLastItem[0] + 1;

    //cerr << "peer 0:" << sum << "\n";

    for(int i = 1; i < params->peers; i++) {
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
    if(sum != params->ni) {
        cerr << "ERROR: ni = "<< params->ni << "!= sum = " << sum << endl;
        printError(DATASET_DIVISION_ERROR, __FUNCTION__);
        return DATASET_DIVISION_ERROR;
    }

    cout.flush();

    return SUCCESS;
}

Params* init() {

    // Initialize the sketch based on user-supplied parameters
    Params *params = nullptr;
    int rows = 0;

    params = new (nothrow) Params;
    if(!params){
        fprintf(stdout,"Memory allocation of sketch data structure failed\n");
        return nullptr;
    }

    params->name_file = DEFAULT_FILENAME;
    int returnValue = getDatasetSize(params->name_file, rows);
    if ( returnValue < 0) {
        return nullptr;
    }
    params->ni = rows;
    params->domainSize = DEFAULT_DOMAIN_SIZE;
    params->graphType = DEFAULT_GRAPH_TYPE;
    params->peers = DEFAULT_PEERS;
    params->fanOut = DEFAULT_FAN_OUT;
    params->convThreshold = DEFAULT_CONVERGENCE_THRESHOLD;
    params->convLimit = DEFAULT_CONVERGENCE_LIMIT;
    params->roundsToExecute = DEFAULT_ROUND_TO_EXECUTE;
    params->outputOnFile = !params->outputFilename.empty();
    params->autoseed = DEFAULT_AUTO_SEED;

    params->offset = DEFAULT_OFFSET;
    params->bin_limit = DEFAULT_BIN_LIMIT;
    params->alpha = DEFAULT_ALPHA;

    return params;
}

int parse(int argc, char **argv, Params* params) {

    /*** parse command-line parameters ***/
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-di") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing domain size parameter." << endl;
                return -2;
            }
            params->domainSize = stol(argv[i]);
        } else if (strcmp(argv[i], "-p") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of peers parameter." << endl;
                return -2;
            }
            params->peers = stoi(argv[i]);
        } else if (strcmp(argv[i], "-f") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing fan-out parameter." << endl;
                return -2;
            }
            params->fanOut = stoi(argv[i]);
        } else if (strcmp(argv[i], "-s") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing seed parameter" << endl;
                return -2;
            }
            params->graphType = stoi(argv[i]);
        } else if (strcmp(argv[i], "-d") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing graph type parameter" << endl;
                return -2;
            }
            params->graphType = stoi(argv[i]);
        } else if (strcmp(argv[i], "-ct") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing convergence tolerance parameter." << endl;
                return -2;
            }
            params->convThreshold = stod(argv[i]);
        } else if (strcmp(argv[i], "-cl") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing # of consecutive rounds in which convergence is satisfied parameter." << endl;
                return -2;
            }
            params->convLimit = stol(argv[i]);
        } else if (strcmp(argv[i], "-out") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing filename for simulation output." << endl;
                return -2;
            }
            params->outputFilename = string(argv[i]);
        } else if (strcmp(argv[i], "-r") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of rounds to execute.\n";
                return -2;
            }
            params->roundsToExecute = stoi(argv[i]);
        } else if (strcmp(argv[i], "-alp") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of accuracy for DDSketch.\n";
                return -2;
            }
            params->alpha = stod(argv[i]);
        } else if (strcmp(argv[i], "-ofs") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of offset for DDSketch.\n";
                return -2;
            }
            params->offset = stoi(argv[i]);
        } else if (strcmp(argv[i], "-bl") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing number of bins limit for DDSketch.\n";
                return -2;
            }
            params->bin_limit = stoi(argv[i]);
        } else if (strcmp(argv[i], "-dt") == 0) {
            i++;
            if (i >= argc) {
                cerr << "Missing file name of dataset.\n";
                return -2;
            }
            params->name_file = argv[i];
        } else if (strcmp(argv[i], "-as") == 0) {
            params->autoseed = true;
        } else {
            usage(argv[0]);
            printError(USAGE_ERROR, __FUNCTION__);
            return USAGE_ERROR;
        }
    }

    return SUCCESS;
}

void usage(char* cmd)
{
    cerr
            << "Usage: " << cmd << "\n"
            << "-def        the default parameters will be loaded\n"
            << "or\n"
            << "Usage: " << cmd << "\n"
            << "-di         domain size\n"
            << "-p          number of peers\n"
            << "-f          fan-out of peers\n"
            << "-s          seed\n"
            << "-d          graph type: 1 geometric 2 Barabasi-Albert 3 Erdos-Renyi 4 regular\n"
            << "-ct         convergence tolerance\n"
            << "-cl         number of consecutive rounds in which convergence must be satisfied\n"
            << "-out        output filename, if specified a file with this name containing all of the peers stats is written\n"
            << "-r          number of rounds to execute\n"
            << "-alp        accuracy for DDSketch\n"
            << "-ofs        offset for DDSketch\n"
            << "-bl         bins limit for DDSketch\n"
            << "-ds         name of dataset\n"
            << "-as         enable autoseeding\n\n";
}

void startTheClock(){
    t1 = high_resolution_clock::now();
}

double stopTheClock() {
    t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    return time_span.count();
}

int getDatasetSize(const string &name_file, int &rows) {

    rows = 0;

    ifstream inputFile(name_file);
    string line;
    if(inputFile.is_open()){
        while (getline(inputFile, line))
            rows++;
    } else {
        printError(FILE_ERROR, __FUNCTION__);
        return FILE_ERROR;
    }

    return SUCCESS;
}

int loadDataset(const string &name_file, double *dataset) {

    if (!dataset) {
        printError(NULL_POINTER_ERROR, __FUNCTION__);
        return NULL_POINTER_ERROR;
    }

    ifstream inputFile(name_file);
    string line;
    int row = 0;

    if(inputFile.is_open()){
        while (getline(inputFile, line)) {

            try {
                dataset[row] = stod(line);
            }
            catch (int e) {
                printError(FILE_ERROR, __FUNCTION__);
                inputFile.close();
                return FILE_ERROR;
            }

            row++;
        }
    } else {
        printError(FILE_ERROR, __FUNCTION__);
        return FILE_ERROR;
    }

    inputFile.close();

    return SUCCESS;
}
int testQuantile(DDS_type *dds, double* stream, long n_element, stringstream &output) {

    output << "\nTest quantile with alpha=" << dds->alpha << endl << endl;

    if (!stream) {
        printError(NULL_POINTER_ERROR, __FUNCTION__);
        return NULL_POINTER_ERROR;
    }

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    int returnValue = -1;

    // Determine the quantile
    float q[] = { 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99};
    int n_q = 11;

    output  << string(60, '-') << endl
            << "|" << setw(10) << "quantile" << "|" << setw(15) << "estimate" << "|" << setw(15) << "real" << "|"  << setw(15)<< "error" << "|" << endl
            << string(60, '-') << endl;

    for ( int i = 0; i < n_q; i++ ) {

        int idx = floor(1+q[i]*double(n_element-1));
        // determine the correct answer
        // i.e., the number stored at the index (idx-1) in the sorted permutation of the vector
        // note that we are not sorting the vector, we are using the quickselect() algorithm
        // which in C++ is available as std::nth_element
        nth_element(stream, stream + (idx-1), stream +  n_element);
        double quantile;
        returnValue = DDS_GetQuantile(dds, q[i], quantile);
        if (returnValue < 0 ) {
            return returnValue;
        }

        double error = abs((quantile-stream[idx-1])/stream[idx-1]);

        output << "|" << setw(10) << q[i]<< "|" << setw(15) << quantile << "|" << setw(15) << stream[idx-1] << "|"  << setw(15)<< error << "|" << endl;

    }

    output << string(60, '-') << endl;

    return  returnValue;
}

string getParameters(Params* params){

    stringstream parameters;
    parameters << "dataset = "    << params->name_file << "\n"
               << "n° points = "  << params->ni << "\n"
               << "graph type = " << printGraphType(params->graphType)
               << "peers = " << params->peers << "\n"
               << "fan-out = " << params->fanOut << "\n"
               << "local convergence tolerance = "<< params->convThreshold << "\n"
               << "number of consecutive rounds in which a peer must locally converge = "<< params->convLimit << "\n"
               << "alpha = " << params->alpha << "\n"
               << "offset = " << params->offset << "\n"
               << "bin_limit = " << params->bin_limit << "\n";

    return parameters.str();
}