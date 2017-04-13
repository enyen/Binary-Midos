#ifndef MIDOS_H
#define MIDOS_H

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <ctime>
using namespace std;


/**
 * @brief The database struct : sturctured data format for input database
 */
struct database{
    bool target; // target value
    int m; // # of entry
    int n; // # of attribute
    int k; // # of sub-group requested
    int F; // quality function

    vector<vector<bool> > dataAttributes; // all attributes of all entry in database
    vector<bool> dataLabel; // label of all entry

    database() : target(true), m(0), n(0), k(0), F(0) {}
};


/**
 * @brief The hypo struct : single hypothesis stump
 */
struct hypo{
    int attribute; // stump attribute
    bool value; // stump value
    bool operator==(hypo &hypo_) const{
        return ((hypo_.attribute==this->attribute) && (hypo_.value==this->value));
    }
    hypo() : attribute(0), value(0) {}
};


/**
 * @brief The group struct : group of hypothesis stump
 */
struct group{
    vector<hypo> groupHypo; // group of stump
    float p; // ext(groupHypo, target) / ext(groupHypo)
    float q; // group quality value
    float g; // ext(groupHypo) / # of entry
    bool operator==(group &group_) const{
        if(this->groupHypo.size()!=group_.groupHypo.size()) return false;
        for(int i=0; i<this->groupHypo.size(); i++)
            if(!(this->groupHypo[i]==group_.groupHypo[i])) return false;
        return true;
    }
    group() : p(0), q(0), g(0) {}
};


/**
 * @brief The midos class : algorithm class
 */
class midos
{
private:
    bool _target;
    int _m;
    int _n;
    int _k;
    int _F;
    float _p0;
    database *_db;

public:
    midos(): _target(true), _m(0), _n(0), _k(0), _F(0), _p0(0), _db(NULL) {}
    midos(database *db);

    void subGroupDiscovery(vector<group> &solution);
    void refinementOperator(group &hypo0, vector<group> &hypos_);
    void add2Solution(vector<group> &solution, group &toAdd, float &minQ, int &minQ_idx);
    void try2Add2Solution(vector<group> &solution, group &toAdd, float &minQ, int &minQ_idx);
    bool pruneCheck(group &hypo0, float minQ);
    void properties(group &hypo0, bool target);
    const float get_p0() {return _p0;}
};

#endif // MIDOS_H
