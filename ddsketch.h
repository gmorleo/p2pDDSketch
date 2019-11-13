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
#include <math.h>
#include <limits>
#include <map>
#include <algorithm>
#include <random>
#include <fstream>

using namespace std;

typedef struct DDS_type{
    /// used to allow the skecth storing both positive, 0 and negative values
    int offset;
    /// maximum number of bins
    int bin_limit;
    /// this parameter defines alpha-accuracy of a q-quantile
    float alpha;
    /// this parameter is defined as (1 + alpha)/(1-alpha)
    float gamma;
    /// this is not a required parameter; it is defined as log(gamma)
    float ln_gamma;
    /// this map implements the bins of DDSketch
    map<int, double> *bins;
    /// this parameter keeps track of the number of items added to the sketch
    double n;
} DDS_type;

typedef  struct DDS_interval {
    /// left limit of the interval
    double min;
    /// right limit of the interval
    double max;
    /// lenght limit of the interval
    double lenght;
} DDS_interval;

typedef struct DDS_split_interval {
    /// percent of elements in the precedent bin
    double precedent;
    /// percent of elements in the next bin
    double next;
} DDS_split_interval;

/**
 * @brief               DDS costructor
 * @param offset        Used to allow the skecth storing both positive, 0 and negative values
 * @param bin_limit     The maximum number of bins
 * @param alpha         The alpha-accuraxy level of q-quantile
 * @return              User-supplied parameters
 */
extern DDS_type *DDS_Init(int offset, int bin_limit, float alpha);

/**
 * @brief               DDS destructor
 * @param dds           User-supplied parameters
 */
extern void DDS_Destroy(DDS_type *dds);
/**
 * @brief               Return the number of bins currently in the sketch
 * @param dds           User-supplied parameters
 * @return              The bins size, (number of bins)
 */
extern long DDS_Size(DDS_type *dds);

/**
 * @brief               Given a value x, getKey returns the bucket index
 * @param dds           User-supplied parameters
 * @param item          The the input value
 * @return              The index of the bucket containing the value x
 */
extern int DDS_GetKey(DDS_type *dds, double item);

/**
 * @brief               Given a value x, getKey returns the bucket ind
 * @param dds           User-supplied parameters
 * @param item          The the input value
 * @param ln_gamma      this is not a required parameter; it is defined as log(gamma)
 * @return              The index of the bucket containing the value x
 */
extern int DDS_GetKey(DDS_type *dds, double item, float ln_gamma);

/**
 * @brief               Given a bucket index (i), this function returns the estimation of the rank x_q
 * @param dds           User-supplied parameters
 * @param i             The key of the bin
 * @return              The estimate of the rank x_q
 */
extern double DDS_GetRank(DDS_type *dds, int i);

/**
 * @brief               This function creates a new bucket with index associated with the value (item), or if that bucket already exists, it simply add 1 to the bucket's counter
 * @param dds           User-supplied parameters
 * @param item          The the input value
 * @return              0 success, -1 error
 */
extern int DDS_Add(DDS_type *dds, double item);

/**
 * @brief               This function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1 otherwise it simply decrements by 1 the bucket's counter
 * @param dds           User-supplied parameters
 * @param item          The the input value
 * @return              0 success
 */
extern int DDS_Delete(DDS_type *dds, double item);

/**
 * @brief               In order to reduce the bucket's number, we need to increase the range of the bucket's index.
 * @param dds           User-supplied parameters
 * @return              0 success, -1 error
 */
extern int DDS_expand(DDS_type *dds);

/**
 * @brief               The function computes the estimate of the desired q-quantile (q)
 * @param dds           User-supplied parameters
 * @param q             The desired q-quantile
 * @return              The estimate of the desired q-quantile
 */
extern double DDS_GetQuantile(DDS_type *dds, float q);

/**
 * @brief               Merge function merges the bins in dds1 with the bins of dds2 dds1 is the result of the merge operation
 * @param dds1          Parameters of the first sketch
 * @param dds2          Parameters of the second sketch
 */
extern void DDS_merge(DDS_type *dds1, DDS_type *dds2);

/**
 * @brief                   This function computes the sum of all counter of the bins
 * @param dds               Parameters of the sketch
 * @return                  int
 */
extern int DDS_SumBins(DDS_type *dds);

/**
 * @brief                   This function prints the bins map in a CSV file
 * @param name              File name
 * @param bins              Bins map
 * @return                  0 success
 */
extern int DDS_PrintCSV(string name, map<int,double> *bins);

/**
 * @brief                   This function checks if all the elements in the stream have a corresponding bucket
 * @param dds               Parameters of the sketch
 * @param item              Input value
 * @return                  0 success, -1 failed
 */
extern int DDS_CheckAll(DDS_type *dds, double item);

/**
 * @brief                   This function computes the range of a key
 * @param dds               Parameters of the sketch
 * @param i                 Key
 * @param gamma             Gamma value according alpha
 * @return                  DDS_interval
 */
extern DDS_interval *DDS_getInterval(DDS_type *dds, int i, float gamma);

/**
 * @brief                   This function computes the percent of elements that are going to be distributed in the precedent or next bin.
 * @param dds               Parameters of the sketch
 * @param i                 Index of old bin of value (x) according the old alpha
 * @param k                 Index of new bin of value (x) according the new alpha
 * @param gamma_i           Gamma value according the old alpha
 * @param gamma_k           Gamma value according the new alpha
 * @return                  DDS_split_interval
 */
extern DDS_split_interval *DDS_getSplitInterval(DDS_type *dds, int i, int k, float gamma_i, float gamma_k);

/**
 * @brief                   This function expands all the bins in the map, increasing alpha by 0.01. The values in the old range are redistributed in the new range with a proportional way (supposing all values uniforming distributed on the interval)
 * @param dds               Parameters of the sketch
 * @return                  0 success
 */
extern int DDS_expandProportional(DDS_type *dds);

/**
 * @brief                   This function finalizes the average consensus algorithm, dividing the mean with the weight
 * @param dds               Parameters of the sketch
 * @param weight            Weight
 * @return                  0 success
 */
extern int DDS_finalizeMerge(DDS_type *dds, double weight);

/**
 * @brief                   This function return true or false if the two sketches are mergeable or not
 * @param dds1              Parameters of first sketch
 * @param dds2              Parameters of second sketch
 * @return                  true/false
 */
extern bool DDS_isMergeable(DDS_type *dds1, DDS_type *dds2);

/**
 * @brief                   This function replace the bins map and parameters ( alpha, gamma, ln_gamma, n ) of the first sketch with the ones in the second sketch
 * @param dds1              Parameters of first sketch
 * @param dds2              Parameters of second sketch
 * @return                  0 success
 */
extern int DDS_replaceBinMap(DDS_type *dds1, DDS_type *dds2);

extern int DDS_prova(DDS_type *dds);
extern int DDS_prova2(DDS_type *dds);