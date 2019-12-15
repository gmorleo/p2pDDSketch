//
// Created by giuseppe on 25/11/19.
//

#include <iostream>

#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define RESET   "\033[0m"

const int SUCCESS = 0;
const int GENERIC_ERROR = -1;
const int MEMORY_ERROR = -2;
const int FILE_ERROR = -3;
const int SKETCH_ERROR = -4;
const int MERGE_ERROR = -5;
const int QUANTILE_ERROR = -6;
const int UNKNOWN_COLLAPSE_TYPE = -7;
const int COPY_ERROR = -8;
const int NULL_POINTER_ERROR = -9;
const int GRAPH_GENERATION_ERROR = -10;
const int UNKNOWN_GRAPH_TYPE = -11;
const int PARAM_ERROR = -12;
const int USAGE_ERROR = -13;
const int CONFLICTING_OPTIONS = -14;
const int DATASET_DIVISION_ERROR = -15;
const int EXIT = -16;

using namespace std;

extern int printError(int error, const string& nameFunction);
