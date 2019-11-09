/********************************************************************
 DDSketch

 An algorithm for tracking quantiles in data streams

 Charles Masson, Jee E. Rim, and Homin K. Lee. 2019. DDSketch: a fast and fully-mergeable quantile sketch with relative-error guarantees. Proc. VLDB Endow. 12, 12 (August 2019), 2195-2205. DOI: https://doi.org/10.14778/3352063.3352135

 This implementation by
 by Giuseppe Morleo
 University of Salento, Italy

 *********************************************************************/
/// \file

#include "ddsketch.h"



DDS_type *DDS_Init(int offset, int bin_limit, float alpha)
{
    
    // Initialize the sketch based on user-supplied parameters
    DDS_type *dds = NULL;
    
    dds = new(DDS_type); // do not use the C malloc() function: it does not call C++ constructors
    if(!dds){
        fprintf(stdout,"Memory allocation of sketch data structure failed\n");
        return NULL;
    }
    
    dds->offset = offset;
    dds->bin_limit = bin_limit;
    dds->alpha = alpha;
    dds->gamma = (1 + dds->alpha)/(1 - dds->alpha);
    dds->ln_gamma = log(dds->gamma);
    dds->bins = new map<int, double>();
    if(!dds->bins){
        fprintf(stdout,"Memory allocation of sketch map failed\n");
        delete dds;
        return NULL;
    }
    
    dds->n = 0;
    
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

long DDS_Size(DDS_type *dds)
{
    // return the number of bins currently in the sketch
    return dds->bins->size();
}

int DDS_GetKey(DDS_type *dds, double item)
{
    // Given a value (item), returns the correspondning bucket index
    int key;
    
    if (item > 0) {
        key = int(ceil((log(item))/dds->ln_gamma)) + dds->offset;
    } else if (item < 0) {
        key = -int(ceil((log(-item))/dds->ln_gamma)) - dds->offset;
    } else {
        key = 0;
    }
    
    return key;
    
}

int DDS_GetKey(DDS_type *dds, double item, float ln_gamma)
{
    // Given a value (item), returns the correspondning bucket index
    int key;
    
    if (item > 0) {
        key = int(ceil((log(item))/ln_gamma)) + dds->offset;
    } else if (item < 0) {
        key = -int(ceil((log(-item))/ln_gamma)) - dds->offset;
    } else {
        key = 0;
    }
    
    return key;
    
}

double DDS_GetRank(DDS_type *dds, int i)
{
    
    // Given a bucket index (i), this function returns the estimation of the rank x_q
    
    if (i > 0) {
        i -= dds->offset;
        return (2 * pow(dds->gamma, i))/(dds->gamma + 1);
    } else {
        i += dds->offset;
        return -(2 * pow(dds->gamma, -i))/(dds->gamma + 1);
    }
}

int DDS_Add(DDS_type *dds, double item)
{
    
    // this function creates a new bucket with index associated with the value (item), or if that bucket already exists, it
    // simply add 1 to the bucket's counter
    
    int key = DDS_GetKey(dds, item);
    (*(dds->bins))[key] += 1;
    dds->n += 1;
    
    if (DDS_Size(dds) > dds->bin_limit ){
        // If the bin size is greater than bin_limit, we need to increase alpha and adapt all of the existing buckets to the new alpha value
       
        int return_value = DDS_expandProportional(dds);
        if(return_value < 0){
            return -1;
        }
            
    }
    
    return 0;
    
}

int DDS_Delete(DDS_type *dds, double item)
{

    // this function deletes the bucket with index associated with the value (item) if it exists and its value is equal to 1
    // otherwise it simply decrements by 1 the bucket's counter

    int key = DDS_GetKey(dds, item);
    map<int, double>::iterator it = dds->bins->find(key);
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
        //dds->n -= 1;
        //cout << "There is no bin associated to item " << item << " with key " << key << endl;

    }

    return 0;

}

int DDS_expand(DDS_type *dds)
{
    // In order to reduce the bucket's number, we need to increase the range of the bucket's index.
    // We compute the new values of gamma and ln_gamma according the new alpha.
    dds->alpha += 0.01;
    cout << "New alpha = " << dds->alpha << endl;
    float new_gamma = (1 + dds->alpha)/(1 - dds->alpha);
    float new_ln_gamma = log(new_gamma);

    double item;
    int key;
    
    // Create new bins map
    map<int,double> *new_bins = NULL;
    new_bins = new map<int, double>();
    if(!new_bins){
        fprintf(stdout,"Memory allocation of a new sketch map failed\n");
        return -1;
    }

    for (auto & bin : (*(dds->bins))) {
        item = DDS_GetRank(dds, bin.first);
        key = DDS_GetKey(dds, item, new_ln_gamma);
        (*new_bins)[key] += bin.second;
    }
    
    // Replace old bins map with new bins map
    dds->bins->swap(*new_bins);
    new_bins->clear();
    delete new_bins;
    
    dds->gamma = new_gamma;
    dds->ln_gamma = new_ln_gamma;
    
    return 0;

}

int DDS_SumBins(DDS_type *dds) {

    int sum = 0;

    for (auto & bin : (*(dds->bins))) {
        sum += bin.second;
    }

    return sum;
}

int DDS_PrintCSV(string name, map<int,double> *bins) {

    ofstream file;
    file.open(name);

    for (auto& b: (*bins)) {
        file << b.first << ", ";
        file << b.second << ",\n";
    }

    file.close();

    return 0;
}


int DDS_CheckAll(DDS_type *dds, double item) {

    int key = DDS_GetKey(dds, item);
    bool check = true;

    auto it = dds->bins->find(key);
    if (it == dds->bins->end()){
        cout << "Not found key = " << key << endl;
        DDS_getInterval(dds, key, dds->gamma);
        check = false;
    }

    return check;
}

double DDS_GetQuantile(DDS_type *dds, float q)
{
    // If the value (q) is not in the [0,1] interval return NaN
    if (q < 0 || q > 1.01) {
        return numeric_limits<double>::quiet_NaN();
    }

    // We need to sum up the buckets until we find the bucket containing the q-quantile value x_q
    auto it = dds->bins->begin();
    int i = it->first;
    int count = it->second;
    
    while (count <= q*(dds->n - 1)) {
        ++it;
        i = it->first;
        count += it->second;
    }
    
    // Return the estimation x_q of bucket index i
    return DDS_GetRank(dds, i);
    
}

void DDS_merge(DDS_type *dds1, DDS_type *dds2)
{

    // Merge function merges the bins in dds1 with the bins of dds2
    // dds1 is the result of the merge operation

    // Check if the bins have the same alpha
    while (!DDS_isMergeable(dds1,dds2)){
        if (dds1->alpha < dds2->alpha) {
            DDS_expandProportional(dds1);
        } else {
            DDS_expandProportional(dds2);
        }
    }

    // Mean
    for ( auto & bin  : (*dds1->bins)) {
        bin.second = bin.second/2;
    }

    // Mean
    for ( auto & bin  : (*dds2->bins)) {
        (*(dds1->bins))[bin.first] += bin.second/2;
    }

    // Check if the new bin size is greater than bin limit
    if ( dds1->bins->size() > dds1->bin_limit ) {
        DDS_expandProportional(dds1);
    }

    // Mean of n
    dds1->n = (dds1->n/2 + dds2->n/2);
    //cout << "Size after merge = " << DDS_Size(dds1) << endl;
}

bool DDS_isMergeable(DDS_type *dds1, DDS_type *dds2) {

    return abs(dds1->alpha - dds2->alpha) <= 0.005;

}

DDS_interval *DDS_getInterval(DDS_type *dds, int i, float gamma) {

    DDS_interval *interval = NULL;
    interval = new DDS_interval{0, 0};

    if ( i > 0) {
        i -= dds->offset;
        interval->min = pow(gamma,i-1);
        interval->max = pow(gamma,i);
    } else {
        i += dds->offset;
        interval->min = -pow(gamma,-(i-1));
        interval->max = -pow(gamma,-i);
    }

    interval->lenght = abs(interval->max-interval->min);

    return interval;
}

DDS_split_interval *DDS_getSplitInterval(DDS_type *dds, int i, int k, float gamma_i, float gamma_k) {

    DDS_split_interval *split = NULL;
    split = new DDS_split_interval{0, 0};

    DDS_interval *old_interval = DDS_getInterval(dds, i, gamma_i);
    DDS_interval *new_interval = DDS_getInterval(dds, k, gamma_k);

    //cout << "old " << old_interval->min << " - " << old_interval->max << endl;
    //cout << "new " << new_interval->min << " - " << new_interval->max << endl;

    // If the old max is greater than the new max, a portion of elements must be distributed in the next bin
    if ( old_interval->max > new_interval->max) {
        split->next = abs((old_interval->max-new_interval->max)/old_interval->lenght);
    }

    // If the old min is fewer than the new min, a portion of elements must be distributed in the precedent bin
    if ( old_interval->min < new_interval->min ) {
        split->precedent = abs((new_interval->min-old_interval->min)/old_interval->lenght);
    }

    // Saturation
    if ( split->next > 1 ) {
        split->next = 1;
    }

    // Saturation
    if ( split->precedent > 1 ) {
        split->precedent = 1;
    }

    //cout << "prec " << split->precedent << " next " << split->next << endl;

    delete new_interval;
    delete old_interval;

    return  split;
}

int DDS_expandProportional(DDS_type *dds){

    // In order to reduce the bucket's number, we need to increase the range of the bucket's index.
    // We compute the new values of gamma and ln_gamma according the new alpha.
    dds->alpha += 0.01;
    float new_gamma = ((1 + dds->alpha)/(1 - dds->alpha));
    float new_ln_gamma = log(new_gamma);
    //cout << "New alpha = " << dds->alpha << endl;

    DDS_split_interval *split = NULL;

    double item;
    int temp;
    int key;

    // Create new bins map
    map<int,double> *new_bins = NULL;
    new_bins = new map<int, double>();
    if(!new_bins){
        fprintf(stdout,"Memory allocation of a new sketch map failed\n");
        return -1;
    }

    for (auto & bin : (*dds->bins)) {

        item = DDS_GetRank(dds, bin.first);
        key = DDS_GetKey(dds, item, new_ln_gamma);
        split = DDS_getSplitInterval(dds, bin.first, key, dds->gamma, new_gamma);

        // Element in precedent bin
        temp = floor(bin.second * split->precedent);
        if ( temp != 0 ) {
            (*new_bins)[key - 1] += temp;
        }

        // Elementi in next bin
        temp = floor(bin.second * split->next);
        if ( temp != 0 ) {
            (*new_bins)[key + 1] += temp;
        }

        // Elementi in current bin
        temp = bin.second - floor(bin.second * split->precedent) - floor(bin.second * split->next);
        if ( temp != 0 ) {
            (*new_bins)[key] += temp;
        }

        delete split;
    }

    // Replace old bins map with new bins map
    dds->bins->swap(*new_bins);
    new_bins->clear();
    delete new_bins;

    //cout << "Size after expand = " << DDS_Size(dds) << endl;

    dds->gamma = new_gamma;
    dds->ln_gamma = new_ln_gamma;

    return 0;
}

int DDS_finalizeGossip(DDS_type *dds, double weight) {

    for(auto & bin: (*dds->bins)) {
        bin.second = ceil(bin.second/weight);
    }
    dds->n = dds->n/weight;

    return 0;
}

int DDS_replaceBinMap(DDS_type *dds1, DDS_type *dds2) {

    dds2->alpha = dds1->alpha;
    dds2->gamma = dds1->gamma;
    dds2->ln_gamma = dds2->ln_gamma;
    dds2->n = dds2->n;

    dds2->bins->clear();
    dds2->bins->insert(dds1->bins->begin(),dds1->bins->end());

    return 0;
}