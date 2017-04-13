#include "midos.h"

/**
 * @brief midos::midos : initialise parameters of algorithm
 * @param db (in) : database passed to this algorithm
 */
midos::midos(database *db){
    _target = db->target;
    _m = db->m;
    _n = db->n;
    _k = db->k;
    _F = db->F;
    int sum = 0;
    for(int i=0; i<db->dataLabel.size(); i++)
        sum += db->dataLabel[i];
    _p0 = _target ? (sum * 1.f / _m) : (1-(sum * 1.f / _m));
    _db = db;
}


/**
 * @brief midos::subGroupDiscovery : main function for discovering sub group
 * @param solution (out) : collection of best k hypothesis
 */
void midos::subGroupDiscovery(vector<group> &solution){
    cout<<"Discovering..."<<endl;
    solution.clear();
    vector<group> groupPools; // working pool

    // initialise first tree layer
    group emptyGroup;
    float minQ = 100;
    int minQ_idx = 0;
    refinementOperator(emptyGroup, groupPools);

    // subsequence tree layer loop
    while(groupPools.size()!=0){
        cout<<"new layer size "<<groupPools.size()<<endl;
        vector<group> newGroupPool; // next layer working pool
        // each leaf in this tree layer loop - breath first search
        for(int i=0; i<groupPools.size(); i++){
            if(solution.size() < _k) // solution pool not yet full
                add2Solution(solution, groupPools[i], minQ, minQ_idx);
            else // solution pool full
                try2Add2Solution(solution, groupPools[i], minQ, minQ_idx);

            if(!pruneCheck(groupPools[i], minQ)){ // check expansion worthiness of leaf
                vector<group> thisGroupPool;
                refinementOperator(groupPools[i], thisGroupPool);
                for(auto thisgroup : thisGroupPool) newGroupPool.push_back(thisgroup); // collect new leaf
            }
        }
        groupPools.clear();
        for(auto thisgroup : newGroupPool) groupPools.push_back(thisgroup); // new leaf pool to working pool
    }
    cout<<"Sub-groups discovered."<<endl<<endl;
}


/**
 * @brief midos::refinementOperator : generate child based on parent hypothesis
 * @param hypo0 (in) : parent hypothesis
 * @param hypos_ (out) : child hypothesis
 */
void midos::refinementOperator(group &hypo0, vector<group> &hypos_){
    hypos_.clear();
    group newgroup;
    for(auto thishypo : hypo0.groupHypo) newgroup.groupHypo.push_back(thishypo); // initialise working pool with parent hypo
    int ord0 = (hypo0.groupHypo.size()) ? hypo0.groupHypo[hypo0.groupHypo.size()-1].attribute : -1; // search for highest attribute order
    // subsequence attribute loop
    for(int i=ord0+1; i<_n; i++){
        hypo newhypo;
        newhypo.attribute = i;
        // position new hypo
        newhypo.value = true;
        newgroup.groupHypo.push_back(newhypo); // add to working working pool
        properties(newgroup, _target); // compute properties
        hypos_.push_back(newgroup); // add to output pool
        newgroup.groupHypo.pop_back(); // remove new hypo from working pool
        // negative new hypo
        newhypo.value = false;
        newgroup.groupHypo.push_back(newhypo);
        properties(newgroup, _target);
        hypos_.push_back(newgroup);
        newgroup.groupHypo.pop_back();
    }
}


/**
 * @brief midos::add2Solution : add hypothesis to solution pool when there is still space
 * @param solution (out) : solution pool
 * @param toAdd (in) : new hypothesis to be added
 * @param minQ (in) : current min Qvalue, (out) : new min Qvalue
 * @param minQ_idx (in) : current min Qvalue index, (out) : new min Qvalue index
 */
void midos::add2Solution(vector<group> &solution, group &toAdd, float &minQ, int &minQ_idx){
    // update min Qvalue and index
    minQ = (toAdd.q < minQ) ? toAdd.q : minQ;
    minQ_idx = (toAdd.q < minQ) ? (solution.size()) : minQ_idx;
    solution.push_back(toAdd); // add to solution pool
}


/**
 * @brief midos::try2Add2Solution : if better than the worst, remove the worst and add the new hypo to solution pool
 * @param solution (out) : solution pool
 * @param toAdd (in) : new hypothesis to be tested
 * @param minQ (in) : current min Qvalue, (out) : new min Qvalue
 * @param minQ_idx (in) : current min Qvalue index, (out) : new min Qvalue index
 */
void midos::try2Add2Solution(vector<group> &solution, group &toAdd, float &minQ, int &minQ_idx){
    if(toAdd.q > minQ){
    // better than the worst in solution pool
        solution.erase(solution.begin() + minQ_idx); // remove the worst
        solution.push_back(toAdd); // add new hypo to solution pool
        // update min Qvalue and index
        minQ = solution[0].q;
        minQ_idx = 0;
        for(int i=1; i<solution.size(); i++){
            minQ = (solution[i].q < minQ) ? solution[i].q : minQ;
            minQ_idx = (solution[i].q < minQ) ? i : minQ_idx;
        }
    }
}


/**
 * @brief midos::pruneCheck : check if leaf should be pruned
 * @param hypo0 (in) : hypo to be tested
 * @param minQ (in) : current min Qvalue
 * @return (True) : prune, (False) : no prune
 */
bool midos::pruneCheck(group &hypo0, float minQ){
    float qmax;
    // respective qmax formular
    if(_F == 1)
        qmax = sqrt(hypo0.g) * max(_p0, (1-_p0));
    else if(_F == 2)
        qmax = hypo0.g / (1-hypo0.g) * pow(max(_p0, (1-_p0)), 2);
    else if(_F == 3)
        qmax = hypo0.g + 1 - _p0;
    return (qmax < minQ);
//    return ((qmax < minQ) || (hypo0.g < 0.15));
}


/**
 * @brief midos::properties : compute g, p, q of hypo
 * @param hypo0 (in) : hypo to be tested
 * @param target (in) : target label of algorithm
 */
void midos::properties(group &hypo0, bool target){
    int ext_tar = 0;
    int ext = 0;
    // test every entry in database
    for(int i=0; i<_m; i++){
        int pass = 1;
        // match entry with hypo
        for(int j=0; j<hypo0.groupHypo.size(); j++)
            if(_db->dataAttributes[i][hypo0.groupHypo[j].attribute] != hypo0.groupHypo[j].value){
                pass = 0; // unmatched
                break;
            }
        // update ext(hypo) and ext(hypo, target)
        ext += pass;
        ext_tar += (_db->dataLabel[i] == target) ? pass : 0;
    }
    // compute g, p, and q with respective formular
    hypo0.g = 1.f * ext / _m;
    hypo0.p = 1.f * ext_tar / (ext+0.001);
    if(_F == 1)
        hypo0.q = sqrt(hypo0.g) * max(hypo0.p-_p0, _p0-hypo0.p);
    else if(_F == 2)
        hypo0.q = hypo0.g / (1-hypo0.g) * pow(hypo0.p-_p0, 2);
    else if(_F == 3)
        hypo0.q = hypo0.g*(2*hypo0.p-1) + 1 - _p0;
}
