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
        fprintf(stdout,"Memory allocation of sketch map failed\n");
        delete dds;
        return nullptr;
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

int DDS_NewKey(DDS_type* dds,  double i, int of){

    if (i > 0) {
        i -= dds->offset;
        i = ceil((i+of)/2);
        i += dds->offset;
    } else {
        i += dds->offset;
        i = floor((i+of)/2);
        i -= dds->offset;
    }

    return int(i);
}

double DDS_GetBound(DDS_type *dds, int i) {

    // if the bucket index is positive subtract the offset, otherwise add the offset
    // then, compute the bound by exponentiating gamma

    double bound;

    if ( i > 0) {
        i -= dds->offset;
        bound = pow(dds->gamma,i);
    } else {
        i += dds->offset;
        bound = -pow(dds->gamma,-i);
    }

    return bound;

}

double DDS_GetBound(DDS_type *dds, int i, double gamma) {

    // if the bucket index is positive subtract the offset, otherwise add the offset
    // then, compute the bound by exponentiating gamma

    double bound;

    if ( i > 0) {
        i -= dds->offset;
        bound = pow(gamma,i);
    } else {
        i += dds->offset;
        bound = -pow(gamma,-i);
    }

    return bound;

}

int DDS_AddCollapse(DDS_type *dds, double item)
{

    // this function creates a new bucket with index associated with the value (item), or if that bucket already exists, it
    // simply add 1 to the bucket's counter

    int key = DDS_GetKey(dds, item);
    (*(dds->bins))[key] += 1;
    dds->n += 1;

    if (DDS_Size(dds) > dds->bin_limit ){
        // If the bin size is greater than bin_limit, we need to increase alpha and adapt all of the existing buckets to the new alpha value

        // collapse with gamma^2
        int return_value = DDS_Collapse(dds);
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

int DDS_Collapse(DDS_type *dds) {

    //cout << "Size before collapse: " << DDS_Size(dds) << " sum: " << DDS_SumBins(dds) << endl;

    // determine new gamma and alpha values to be used

    dds->gamma = pow(dds->gamma, 2);
    dds->ln_gamma = log(dds->gamma);
    dds->alpha = (2 * dds->alpha) / (1 + pow(dds->alpha, 2));

    // Create new bins map
    map<int, double> *new_bins = NULL;
    new_bins = new(nothrow) map<int, double>();
    if (!new_bins) {
        fprintf(stdout, "Memory allocation of a new sketch map failed\n");
        return -1;
    }

    // scan the buckets
    for (auto it = dds->bins->begin(); it != dds->bins->end(); ++it) {

        int key = it->first;
        // handle positive bucket indexes
        if ( key > 0 ) {

            // if the bucket index is even
            if ( key%2 == 0 ) {

                // determine the new bucket index corresponding to the collapsing of this bucket and the preceding bucket (whose key is odd)
                int new_key = DDS_NewKey(dds, key, -1);
                (*new_bins)[new_key] += it->second;

            } else {

                // the bucket index is odd
                // determine the new bucket index corresponding to the collapsing of this bucket and the next bucket (whose key is even)
                // check if the next bucket exists
                int new_key = DDS_NewKey(dds, key, +1);
                auto next_bin = next(it, 1);

                if ( next_bin->first == key+1) {
                    (*new_bins)[new_key] += it->second + next_bin->second;
                    // we have to skip the next bucket because we have already considered it in this interaction
                    ++it;
                    if ( it == dds->bins->end() ) {
                        break;
                    }
                } else {
                    (*new_bins)[new_key] += it->second;
                }
            }
        } else {
            // handle negative bucket indexes

            // if the bucket index is even
            if ( key%2 == 0 ) {

                // determine the new bucket index corresponding to the collapsing of this bucket and the next bucket (whose key is odd)
                // check if the next bucket exist
                int new_key = DDS_NewKey(dds, key, +1);
                auto next_bin = next(it, 1);

                if ( next_bin->first == key+1) {
                    (*new_bins)[new_key] += it->second + next_bin->second;
                    // we have to skip the next bucket because we have already considered it in this interaction
                    ++it;
                    if ( it == dds->bins->end() ) {
                        break;
                    }
                } else {
                    (*new_bins)[new_key] += it->second;
                }
            } else {
                int new_key = DDS_NewKey(dds, key, -1);
                (*new_bins)[new_key] += it->second;
            }
        }
    }

    // Replace old bins map with new bins map
    dds->bins->swap(*new_bins);
    new_bins->clear();
    delete new_bins;

    //cout << "Size after collapse = " << DDS_Size(dds) << " sum: " << DDS_SumBins(dds) << " alpha: " << dds->alpha << " gamma: " << dds->gamma << endl << endl;

    return 0;
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
            DDS_Collapse(dds1);
        } else {
            DDS_Collapse(dds2);
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
        DDS_Collapse(dds1);
    }

    // Mean of n
    dds1->n = (dds1->n/2 + dds2->n/2);
    //cout << "Size after merge = " << DDS_Size(dds1) << endl;
}

bool DDS_isMergeable(DDS_type *dds1, DDS_type *dds2) {

    return abs(dds1->alpha - dds2->alpha) <= 0.005;

}

int DDS_finalizeMerge(DDS_type *dds, double weight) {

    for(auto & bin: (*dds->bins)) {
        bin.second = (bin.second/weight);
    }
    dds->n = dds->n/weight;

    return 0;
}

int DDS_replaceBinMap(DDS_type *dds1, DDS_type *dds2) {

    dds2->alpha = dds1->alpha;
    dds2->gamma = dds1->gamma;
    dds2->ln_gamma = dds1->ln_gamma;
    dds2->n = dds1->n;

    dds2->bins->clear();
    dds2->bins->insert(dds1->bins->begin(),dds1->bins->end());

    return 0;
}

int DDS_SumBins(DDS_type *dds) {

    double sum = 0;

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

    auto it = dds->bins->find(key);
    if (it == dds->bins->end()){
        cout << "Not found key = " << key << endl;
    }

    return 0;
}