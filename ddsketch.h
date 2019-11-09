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
    double min;
    double max;
    double lenght;
} DDS_interval;

typedef struct DDS_split_interval {
    double precedent;
    double next;
} DDS_split_interval;

/**
 * /brief               DDS costructor
 * @param offset        Used to allow the skecth storing both positive, 0 and negative values
 * @param bin_limit     The maximum number of bins
 * @param alpha         The alpha-accuraxy level of q-quantile
 * @return              User-supplied parameters
 */
extern DDS_type *DDS_Init(int offset, int bin_limit, float alpha);

/**
 * /brief               DDS destructor
 * @param dds           User-supplied parameters
 */
extern void DDS_Destroy(DDS_type *dds);
/**
 * \brief               Return the number of bins currently in the sketch
 * @param dds           User-supplied parameters
 * @return              The bins size, (number of bins)
 */
extern long DDS_Size(DDS_type *dds);

/**
 * \brief               Given a value x, getKey returns the bucket index
 * @param dds           User-supplied parameters
 * @param item          The the input value
 * @return              The index of the bucket containing the value x
 */
extern int DDS_GetKey(DDS_type *dds, double item);

/**
 * \brief               Given a value x, getKey returns the bucket ind
 * @param dds           User-supplied parameters
 * @param item          The the input value
 * @param ln_gamma      this is not a required parameter; it is defined as log(gamma)
 * @return              The index of the bucket containing the value x
 */
extern int DDS_GetKey(DDS_type *dds, double item, float ln_gamma);

/**
 * \brief               Given a bucket index (i), this function returns the estimation of the rank x_q
 * @param dds           User-supplied parameters
 * @param i             The key of the bin
 * @return              The estimate of the rank x_q
 */
extern double DDS_GetRank(DDS_type *dds, int i);

/**
 * \brief               This function creates a new bucket with index associated with the value (item), or if that bucket already exists, it simply add 1 to the bucket's counter
 * @param dds           User-supplied parameters
 * @param item          The the input value
 * @return              0 success, -1 error
 */
extern int DDS_Add(DDS_type *dds, double item);

/**
 * \brief               This function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1 otherwise it simply decrements by 1 the bucket's counter
 * @param dds           User-supplied parameters
 * @param item          The the input value
 * @return              0 success
 */
extern int DDS_Delete(DDS_type *dds, double item);

/**
 * \brief               In order to reduce the bucket's number, we need to increase the range of the bucket's index.
 * @param dds           User-supplied parameters
 * @return              0 success, -1 error
 */
extern int DDS_expand(DDS_type *dds);

/**
 * \brief               The function computes the estimate of the desired q-quantile (q)
 * @param dds           User-supplied parameters
 * @param q             The desired q-quantile
 * @return              The estimate of the desired q-quantile
 */
extern double DDS_GetQuantile(DDS_type *dds, float q);

/**
 * \brief               Merge function merges the bins in dds1 with the bins of dds2 dds1 is the result of the merge operation
 * @param dds1          User-supplied parameters
 * @param dds2          User-supplied parameters
 */
extern void DDS_merge(DDS_type *dds1, DDS_type *dds2);

/**
 *
 * @param dds
 * @return
 */
extern int DDS_SumBins(DDS_type *dds);

/**
 *
 * @param dds
 * @return
 */
extern int DDS_PrintCSV(DDS_type *dds);

/**
 *
 * @param dds
 * @param item
 * @return
 */
extern int DDS_CheckAll(DDS_type *dds, double item);

/**
 *
 * @param dds
 * @param i
 * @param gamma
 * @return
 */
extern DDS_interval *DDS_getInterval(DDS_type *dds, int i, float gamma);

/**
 *
 * @param dds
 * @param i
 * @param k
 * @param gamma_i
 * @param gamma_k
 * @return
 */
extern DDS_split_interval *DDS_getSplitInterval(DDS_type *dds, int i, int k, float gamma_i, float gamma_k);

/**
 *
 * @param dds
 * @return
 */
extern int DDS_expandProportional(DDS_type *dds);

extern int DDS_expandProportional(DDS_type *dds1, DDS_type *dds2, map<int,double> *bins);

/**
 *
 * @param dds1
 * @param dds2
 * @return
 */
extern int DDS_mergeGossip(DDS_type *dds1, DDS_type *dds2);