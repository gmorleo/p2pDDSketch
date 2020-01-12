/********************************************************************
DDSketch

An algorithm for tracking quantiles in data streams

Charles Masson, Jee E. Rim, and Homin K. Lee. 2019. DDSketch: a fast and fully-mergeable quantile sketch with relative-error guarantees. Proc. VLDB Endow. 12, 12 (August 2019), 2195-2205. DOI: https://doi.org/10.14778/3352063.3352135

This implementation by
by Giuseppe Morleo
University of Salento, Italy

*********************************************************************/
/// \file

#include <fstream>
#include <iomanip>
#include "ddsketch.h"
#include "error.h"

#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define RESET   "\033[0m"

DDS_type *DDS_Init(int offset, int bin_limit, double alpha)
{

    // Initialize the sketch based on user-supplied parameters
    DDS_type *dds = nullptr;

    dds = new (nothrow) (DDS_type); // do not use the C malloc() function: it does not call C++ constructors
    if(!dds){
        fprintf(stdout,"Memory allocation of sketch data structure failed\n");
        return nullptr;
    }

    dds->offset = offset;
    dds->bin_limit = bin_limit;
    dds->alpha = alpha;
    dds->gamma = (1 + dds->alpha)/(1 - dds->alpha);
    dds->ln_gamma = log(dds->gamma);
    dds->bins = new map<int, double>();
    if(!dds->bins){
        printError(MEMORY_ERROR, __FUNCTION__);
        delete dds;
        return nullptr;
    }

    dds->n = 0;

    dds->max = numeric_limits<int>::min();
    dds->min = numeric_limits<int>::max();

    dds->min_value = pow(dds->gamma,pow(2,29));

    return dds;
}

void DDS_Destroy(DDS_type *dds)
{
    // get rid of a sketch and free up the space
    if (!dds) return;

    // deallocate the map
    if (dds->bins){
        delete dds->bins;
    }

    // now free the whole data structure
    delete dds;
}

int DDS_Size(DDS_type *dds, int &size)
{
    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // return the number of bins currently in the sketch
    size = dds->bins->size();

    return SUCCESS;
}

int DDS_GetKey(DDS_type *dds, double item, int &key)
{

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // Given a value (item), returns the correspondning bucket index

    if (item > 0) {
        key = int(ceil((log(item))/dds->ln_gamma)) + dds->offset;
    } else if (item < 0) {
        key = -int(ceil((log(-item))/dds->ln_gamma)) - dds->offset;
    } else if ( abs(item) < dds->min_value) {
        key = 0;
    }

    return SUCCESS;

}

int DDS_GetRank(DDS_type *dds, int key, double &rank) {

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // Given a bucket index (i), this function returns the estimation of the rank x_q

    if (key > 0) {
        key -= dds->offset;
        rank =  (2 * pow(dds->gamma, key))/(dds->gamma + 1);
    } else if (key < 0){
        key += dds->offset;
        rank =  -(2 * pow(dds->gamma, -key))/(dds->gamma + 1);
    } else if ( key == 0) {
        rank = 0;
    }

    return SUCCESS;

}

int DDS_GetValue(DDS_type *dds, int key, double &value) {

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // if the bucket index is positive subtract the offset, otherwise add the offset
    // then, compute the bound by exponentiating gamma
    double bound;

    if ( key > 0) {
        key -= dds->offset;
        value = pow(dds->gamma,key);
    } else {
        key += dds->offset;
        value = -pow(dds->gamma,-key);
    }

    return SUCCESS;

}

int DDS_GetBounds(DDS_type *dds, int k, int i, double &min, double &max) {

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // if the bucket index is positive subtract the offset, otherwise add the offset
    // then, compute the bound by exponentiating gamma

    if ( i > 0) {
        i -= dds->offset;
        max = pow(dds->gamma,i);
    } else {
        i += dds->offset;
        i = -i;
        max = -pow(dds->gamma,i);
    }

    if ( k > 0) {
        k -= dds->offset;
        min = pow(dds->gamma,k);
    } else {
        k += dds->offset;
        k = -k;
        min = -pow(dds->gamma,k);
    }

    return SUCCESS;

}

int DDS_CollapseKey(DDS_type* dds, double i, int of, int &collapse_key){

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // if the bucket index is positive subtract the offset, otherwise add the offset
    // then, compute the new key

    if (i > 0) {
        i -= dds->offset;
        i = ceil((i+of)/2);
        i += dds->offset;
    } else if ( i< 0 ){
        i += dds->offset;
        i = -i;
        i = -ceil((i+of)/2);
        i -= dds->offset;
    } else if ( abs(i) < dds->min_value) {
        i = 0;
    }

    collapse_key = int(i);

    return SUCCESS;
}


int DDS_AddCollapse(DDS_type *dds, double item) {

    int returnValue = -1;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // this function creates a new bucket with index associated with the value (item), or if that bucket already exists, it
    // simply add 1 to the bucket's counter

    int key;
    returnValue = DDS_GetKey(dds, item, key);
    if (returnValue) {
        return returnValue;
    }

    (*(dds->bins))[key] += 1;
    dds->n += 1;

    int size;
    returnValue = DDS_Size(dds, size);
    if (returnValue) {
        return returnValue;
    }

    while ( size > dds->bin_limit ) {
        // While the bin size is greater than binLimit, we need to increase alpha and adapt all of the existing buckets to the new alpha value

        // collapse with gamma^2
        returnValue = DDS_Collapse(dds);
        if(returnValue){
            return returnValue;
        }

        returnValue = DDS_Size(dds, size);
        if (returnValue) {
            return returnValue;
        }
    }

    return returnValue;

}

int DDS_AddCollapseLastBucket(DDS_type *dds, double item) {

    int returnValue = -1;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // this function creates a new bucket with index associated with the value (item), or if that bucket already exists, it
    // simply add 1 to the bucket's counter

    int key;
    returnValue = DDS_GetKey(dds, item, key);
    if (returnValue) {
        return returnValue;
    }

    (*(dds->bins))[key] += 1;
    dds->n += 1;

    int size;
    returnValue = DDS_Size(dds, size);
    if (returnValue) {
        return returnValue;
    }

    if ( size > dds->bin_limit ) {
        // If the bin size is greater than binLimit, we need to increase alpha and adapt all of the existing buckets to the new alpha value

        // collapse the second last bucket into the last bucket
        returnValue = DDS_CollapseLastBucket(dds);
        if(returnValue){
            return returnValue;
        }

    }

    return returnValue;

}

int DDS_AddCollapseFirstBucket(DDS_type *dds, double item) {

    int returnValue = -1;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // this function creates a new bucket with index associated with the value (item), or if that bucket already exists, it
    // simply add 1 to the bucket's counter

    int key;
    returnValue = DDS_GetKey(dds, item, key);
    if (returnValue) {
        return returnValue;
    }

    (*(dds->bins))[key] += 1;
    dds->n += 1;

    int size;
    returnValue = DDS_Size(dds, size);
    if (returnValue) {
        return returnValue;
    }

    if ( size > dds->bin_limit ) {
        // If the bin size is greater than binLimit, we need to increase alpha and adapt all of the existing buckets to the new alpha value

        // collapse the second bucket into the first bucket
        returnValue = DDS_CollapseFirstBucket(dds);
        if(returnValue){
            return returnValue;
        }

    }

    return returnValue;

}

int DDS_DeleteCollapse(DDS_type *dds, double item) {

    int returnValue = -1;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // this function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1
    // otherwise it simply decrements by 1 the bucket's counter

    int key;
    returnValue = DDS_GetKey(dds, item, key);
    if (returnValue) {
        return returnValue;
    }

    auto it = dds->bins->find(key);
    if (it != dds->bins->end()){

        // the bin associate to key actually exists
        // check its value: if it is 1 we erase the bin
        // otherwise we decrement by one the bin

        if(it->second == 1){
            dds->bins->erase(it);
            dds->n -= 1;
            //cout << "Deleted bin associated to item " << item << endl;
        }
        else{
            (*(dds->bins))[key] -= 1;
            dds->n -= 1;
            //cout << "Decremented bin associated to item " << item << endl;
        }


    }
    else{
        returnValue = DDS_RemoveOffset(dds, key);
        if (returnValue) {
            return returnValue;
        }
        //dds->n -= 1;
        //cout << "There is no bin associated to item " << item << " with key " << key << endl;

    }

    return returnValue;

}

int DDS_DeleteCollapseLastBucket(DDS_type *dds, double item) {

    int returnValue = -1;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // this function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1
    // otherwise it simply decrements by 1 the bucket's counter

    int key;
    returnValue = DDS_GetKey(dds, item, key);
    if (returnValue) {
        return returnValue;
    }

    map<int, double>::iterator it;

    if ( key >= dds->min && key <= dds->max ) {
        // the key within the [min,max] interval
        auto last_bucket = dds->bins->rbegin();
        key = last_bucket->first;
        it = dds->bins->find(key);
    } else {
        // the key is outside the [min,max] interval
        it = dds->bins->find(key);
    }

    if (it != dds->bins->end()){

        // the bin associate to key actually exists
        // check its value: if it is 1 we erase the bin
        // otherwise we decrement by one the bin

        if(it->second == 1){
            dds->bins->erase(it);
            dds->n -= 1;
            //cout << "Deleted bin associated to item " << item << endl;
        }
        else{
            (*(dds->bins))[key] -= 1;
            dds->n -= 1;
            //cout << "Decremented bin associated to item " << item << endl;
        }


    }
    else{
        returnValue = DDS_RemoveOffset(dds, key);
        if (returnValue) {
            return returnValue;
        }
        //dds->n -= 1;
        //cout << "There is no bin associated to item " << item << " with key " << key << endl;

    }

    return returnValue;

}

int DDS_DeleteCollapseFirstBucket(DDS_type *dds, double item) {

    int returnValue = -1;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // this function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1
    // otherwise it simply decrements by 1 the bucket's counter

    int key;
    returnValue = DDS_GetKey(dds, item, key);
    if (returnValue) {
        return returnValue;
    }

    map<int, double>::iterator it;

    if ( key >= dds->min && key <= dds->max ) {
        // the key within the [min,max] interval
        it = dds->bins->begin();
        key = it->first;
    } else {
        // the key is outside the [min,max] interval
        it = dds->bins->find(key);
    }

    if (it != dds->bins->end()){

        // the bin associate to key actually exists
        // check its value: if it is 1 we erase the bin
        // otherwise we decrement by one the bin

        if(it->second == 1){
            dds->bins->erase(it);
            dds->n -= 1;
            //cout << "Deleted bin associated to item " << item << endl;
        }
        else{
            (*(dds->bins))[key] -= 1;
            dds->n -= 1;
            //cout << "Decremented bin associated to item " << item << endl;
        }


    }
    else{
        returnValue = DDS_RemoveOffset(dds, key);
        if (returnValue) {
            return returnValue;
        }
        //dds->n -= 1;
        //cout << "There is no bin associated to item " << item << " with key " << key << endl;

    }

    return returnValue;

}

int DDS_GetQuantile(DDS_type *dds, float q, double &quantile) {

    int returnValue = -1;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    if (q < 0 || q > 1.01) {
        printError(QUANTILE_ERROR, __FUNCTION__);
        return QUANTILE_ERROR;
    }

    // We need to sum up the buckets until we find the bucket containing the q-quantile value x_q
    auto it = dds->bins->begin();

    int key = it->first;
    double count = it->second;
    double stop = q*double(dds->n - 1);

    while ( count <= stop) {

        ++it;
        key = it->first;
        count += it->second;

    }

    // Return the estimation x_q of bucket index i
    returnValue = DDS_GetRank(dds, key, quantile);
    if (returnValue) {
        return returnValue;
    }

    return returnValue;

}

int DDS_MergeCollapse(DDS_type *dds1, DDS_type *dds2) {

    int returnValue = -1;

    if(!dds1 || !dds2){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    int size1, size2;

    returnValue = DDS_Size(dds1, size1);
    if (returnValue) {
        return returnValue;
    }

    returnValue = DDS_Size(dds2, size2);
    if (returnValue) {
        return returnValue;
    }
    //cout << "Size before merge sketch1 = " << size1 << " sketch2 = " << size2 << endl;

    // Check if the bins have the same alpha
    while (fabs(dds1->alpha - dds2->alpha) > 0.001){
        if (dds1->alpha < dds2->alpha) {
            //cout << endl << BOLDRED << "Collapsing first sketch..." << RESET << endl;
            returnValue = DDS_Collapse(dds1);
            if (returnValue) {
                return returnValue;
            }
        } else {
            //cout << endl << BOLDRED << "Collapsing second sketch..." << RESET << endl;
            returnValue = DDS_Collapse(dds2);
            if (returnValue) {
                return returnValue;
            }
        }
    }
    
    //cout << endl << BOLDRED << "Size of first sketch before merge: " << dds1->bins->size() << RESET << endl;
    //cout << endl << BOLDRED << "Size of second sketch before merge: " << dds2->bins->size() << RESET << endl;

    // Merge function merges the bins in dds1 with the bins of dds2
    // dds1 is the result of the merge operation
/*    for (auto received_bin : (*(dds2->bins))) {
        (*(dds1->bins))[received_bin.first] += received_bin.second;
        dds1->n += received_bin.second;
    }*/

    for ( auto & bin  : (*dds2->bins)) {
        (*(dds1->bins))[bin.first] += bin.second;
    }

    // Mean
    for ( auto & bin  : (*dds1->bins)) {
        bin.second = bin.second/2;
    }
    
    //cout << endl << BOLDRED << "Size of first sketch after merge but before optional collapse: " << dds1->bins->size() << RESET << endl;

    // Check if the new bin size is greater than bin limit
    returnValue = DDS_Size(dds1, size1);
    if (returnValue < 0) {
        return returnValue;
    }

    while ( size1 > dds1->bin_limit ) {

        // If the bin size is more then the bin limit, we need to collapse using gamma^2 instead of gamma
        //cout << endl << BOLDRED << "Collapsing the merged sketch..." << RESET << endl;
        returnValue = DDS_Collapse(dds1);
        if (returnValue) {
            return returnValue;
        }

        returnValue = DDS_Size(dds1, size1);
        if (returnValue) {
            return returnValue;
        }
    }

    // Mean of n
    dds1->n = (dds1->n/2 + dds2->n/2);
    //cout << "Size after merge = " << DDS_Size(dds1) << endl;

    //cout << "Size after merge = " << size1 << " merge time = " << time << endl << endl;

    return returnValue;
}

int DDS_MergeCollapseLastBucket(DDS_type *dds1, DDS_type *dds2) {

    int returnValue = -1;

    if(!dds1 || !dds2){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    int size1, size2;

    returnValue = DDS_Size(dds1, size1);
    if (returnValue) {
        return returnValue;
    }

    returnValue = DDS_Size(dds2, size2);
    if (returnValue) {
        return returnValue;
    }

    //cout << "Size before merge sketch1 = " << size1 << " sketch2 = " << size2 << endl;

    // Check if the bins have the same alpha
    if (fabs(dds1->alpha - dds2->alpha) > 0.001){
        printError(MERGE_ERROR, __FUNCTION__);
        return MERGE_ERROR;
    }

    // Merge function merges the bins in dds1 with the bins of dds2
    // dds1 is the result of the merge operation
    for (auto received_bin : (*(dds2->bins))) {
        (*(dds1->bins))[received_bin.first] += received_bin.second;
        dds1->n += received_bin.second;
    }

    // Check if the new bin size is greater than bin limit
    returnValue = DDS_Size(dds1, size1);
    if (returnValue) {
        return returnValue;
    }

    if ( dds2->min < dds1->min ) {
        dds1->min = dds2->min;
    }

    if ( dds2->max > dds1->max ) {
        dds1->max = dds2->max;
    }

    while (size1 > dds1->bin_limit){
        // If the bin size is more then the bin limit, we need to collapse the second last bucket in the last bucket

        returnValue = DDS_CollapseLastBucket(dds1);
        if (returnValue) {
            return returnValue;
        }

        returnValue = DDS_Size(dds1, size1);
        if (returnValue) {
            return returnValue;
        }
    }

    //cout << "Size after merge = " << size1 << " merge time = " << time << endl << endl;

    return returnValue;
}

int DDS_MergeCollapseFirstBucket(DDS_type *dds1, DDS_type *dds2) {

    int returnValue = -1;

    if(!dds1 || !dds2){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    int size1, size2;

    returnValue = DDS_Size(dds1, size1);
    if (returnValue) {
        return returnValue;
    }

    returnValue = DDS_Size(dds2, size2);
    if (returnValue) {
        return returnValue;
    }

    //cout << "Size before merge sketch1 = " << size1 << " sketch2 = " << size2 << endl;

    // Check if the bins have the same alpha
    if (fabs(dds1->alpha - dds2->alpha) > 0.001){
        printError(MERGE_ERROR, __FUNCTION__);
        return MERGE_ERROR;
    }

    // Merge function merges the bins in dds1 with the bins of dds2
    // dds1 is the result of the merge operation
    for (auto received_bin : (*(dds2->bins))) {
        (*(dds1->bins))[received_bin.first] += received_bin.second;
        dds1->n += received_bin.second;
    }

    // Check if the new bin size is greater than bin limit
    returnValue = DDS_Size(dds1, size1);
    if (returnValue) {
        return returnValue;
    }

    if ( dds2->min < dds1->min ) {
        dds1->min = dds2->min;
    }

    if ( dds2->max > dds1->max ) {
        dds1->max = dds2->max;
    }

    while (size1 > dds1->bin_limit){

        // If the bin size is more then the bin limit, we need to collapse the second bucket in the first bucket

        returnValue = DDS_CollapseFirstBucket(dds1);
        if (returnValue) {
            return returnValue;
        }

        returnValue = DDS_Size(dds1, size1);
        if (returnValue) {
            return returnValue;
        }
    }

    //cout << "Size after merge = " << size1 << " merge time = " << time << endl << endl;

    return returnValue;
}

int DDS_CollapseLastBucket(DDS_type *dds) {

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // collapse the second bucket into the first bucket
    auto last_bucket = dds->bins->rbegin();
    auto second_last_bucket = ++dds->bins->rbegin();

    if ( second_last_bucket->first < dds->min) {
        dds->min = (second_last_bucket->first);
    }
    if ( last_bucket->first > dds->max) {
        dds->max = last_bucket->first;
    }

    last_bucket->second += second_last_bucket->second;
    dds->bins->erase(second_last_bucket->first);

    return SUCCESS;
}

int DDS_CollapseFirstBucket(DDS_type *dds) {

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    // collapse the second last bucket into the last bucket

    auto first = dds->bins->begin();
    auto second = ++dds->bins->begin();

    if ( first->first < dds->min) {
        dds->min = (first->first);
    }
    if ( second->first > dds->max) {
        dds->max = (second->first);
    }

    first->second += second->second;
    dds->bins->erase(second->first);

    return  SUCCESS;
}

int DDS_Collapse(DDS_type *dds) {

    int returnValue = -1;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    int size;

    returnValue = DDS_Size(dds, size);
    if (returnValue) {
        return returnValue;
    }

    //cout << "Size before collapse: " << size << " alpha: " << dds->alpha << " gamma: " << dds->gamma << endl;

    // determine new gamma and alpha values to be used
    dds->gamma = pow(dds->gamma, 2);
    dds->ln_gamma = log(dds->gamma);
    dds->alpha = (2 * dds->alpha) / (1 + pow(dds->alpha, 2));
    dds->min_value = pow(dds->gamma,pow(2,29));;

    // Create new bins map
    map<int, double> *new_bins = NULL;
    new_bins = new(nothrow) map<int, double>();
    if (!new_bins) {
        printError(MEMORY_ERROR, __FUNCTION__);
        return MEMORY_ERROR;
    }

    // scan the buckets
    for (auto it = dds->bins->begin(); it != dds->bins->end(); ++it) {

        int key = it->first;
        // check if the bucket index is even
        if (key % 2 == 0) {

            int new_key;
            returnValue = DDS_CollapseKey(dds, key, -1, new_key);
            if (returnValue) {
                return returnValue;
            }

            (*new_bins)[new_key] += it->second;

        } else {

            int new_key;
            returnValue = DDS_CollapseKey(dds, key, +1, new_key);
            if (returnValue) {
                return returnValue;
            }

            (*new_bins)[new_key] += it->second;
        }

    }

    // Replace old bins map with new bins map
    dds->bins->swap(*new_bins);

    returnValue = DDS_Size(dds, size);
    if (returnValue) {
        return returnValue;
    }

    //cout << "Size after collapse = " << size << " alpha: " << dds->alpha << " gamma: " << dds->gamma << " collapse time = " << time << endl << endl;

    new_bins->clear();
    delete new_bins;

    return 0;
}

int DDS_PrintCSV(DDS_type* dds, const string& name) {

    int returnValue = -1;
    double min, max;

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return  SKETCH_ERROR;
    }

    ofstream file;
    file.open(name);
    if (file.fail()) {
        printError(FILE_ERROR, __FUNCTION__);
        return FILE_ERROR;
    }

    file << fixed;
    file << setprecision(8);
    file << "key, count, max, min, length, \n";
    for (auto & bin : (*(dds->bins))) {

        returnValue = DDS_GetBounds(dds,bin.first-1, bin.first, min, max);
        if (returnValue) {
            return returnValue;
        }

        int key = bin.first;

        returnValue = DDS_RemoveOffset(dds, key);
        if (returnValue) {
            return returnValue;
        }

        file << key << ", " << bin.second << ", " << max << ", " << min << ", " << max - min << ", \n";

    }

    file.close();

    return returnValue;
}

int DDS_SumBins(DDS_type *dds, int &sum) {

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return  SKETCH_ERROR;
    }

    sum = 0;

    for (auto & bin : (*(dds->bins))) {
        sum += bin.second;
    }

    return SUCCESS;
}

int DDS_RemoveOffset(DDS_type* dds, int &i) {

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    if (i > 0) {
        i -= dds->offset;
    } else {
        i += dds->offset;
    }

    return SUCCESS;
}

int DDS_finalizeGossip(DDS_type *dds, double weight) {

    if(!dds){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    for(auto & bin: (*dds->bins)) {
        bin.second = (bin.second/weight);
    }
    
    dds->n = dds->n/weight;

    return SUCCESS;
}

int DDS_replaceSketch(DDS_type *dds1, DDS_type *dds2) {

    if(!dds1 || !dds2){
        printError(SKETCH_ERROR, __FUNCTION__);
        return SKETCH_ERROR;
    }

    dds2->alpha = dds1->alpha;
    dds2->gamma = dds1->gamma;
    dds2->ln_gamma = dds1->ln_gamma;
    dds2->n = dds1->n;

    dds2->bins->clear();
    dds2->bins->insert(dds1->bins->begin(),dds1->bins->end());

    return SUCCESS;
}