#include "midos.h"

using namespace std;

void loadDatabase(const char *file, database &data);
void saveSubGroup(vector<group> subGroup, float p0, database &db);

int main(int argc, const char *argv[]){
    // check for correct input
    if(argc != 2){
        cout<<"Usage: "<<argv[0]<<" <INPUT>file"<<endl;
        exit(-1);
    }

    database db;
    loadDatabase(argv[1], db);

    midos midos1(&db);
    vector<group> subgroups;

    clock_t tick = clock();
    midos1.subGroupDiscovery(subgroups);
    cout<<"Time spent: "<<(clock()-tick)/(double)CLOCKS_PER_SEC/60.<<" min"<<endl;

    saveSubGroup(subgroups, midos1.get_p0(), db);

    return 0;
}


/**
 * @brief loadDatabase : load database from file
 * @param file (in) : input file name
 * @param data (out) : loaded database struct
 */
void loadDatabase(const char *file, database &data){
    fstream myfile;
    myfile.open(file, ios::in);
    if(!myfile.is_open()){
        cout<<"File not found!"<<endl;
        exit(1);
    }

    // read 4 parameters
    myfile >> data.m;
    myfile >> data.n;
    myfile >> data.k;
    myfile >> data.F;
    char c;
    myfile.get(c);

    // read all entry line
    for(int i=0; i<data.m; i++){
        vector<bool> attributes;
        // read all label and attribute value in the line
        for(int j=0; j<(data.n*2+2); j++){
            myfile.get(c);
            if(((c-'0')==0) || ((c-'0')==1))
                if(j==0) data.dataLabel.push_back((bool)(c-'0'));
                else attributes.push_back(((bool)(c-'0')));
        }
        data.dataAttributes.push_back(attributes);
    }
    myfile.close();
    cout<<endl<<"Database loaded."<<endl;
    cout<<"m="<<data.m<<", n="<<data.n << ", k="<<data.k<<", F="<<data.F<<endl<<endl;

// visualise database
//    for(int i=0; i<data.dataAttributes.size(); i++){
//        cout<<data.dataLabel[i]<<" ";
//        for(int j=0; j<data.dataAttributes[i].size(); j++)
//            cout<<data.dataAttributes[i][j]<< " ";
//        cout<<endl;
//    }
}


/**
 * @brief saveSubGroup : save discovered sub-group to file
 * @param subGroup (in) : all discovered sub-group
 * @param p0 (in) : p0 of the database
 * @param db (in) : original database
 */
void saveSubGroup(vector<group> subGroup, float p0, database &db){
    // sub-group sorting based on Qvalue
    vector<group> subGroupSorted;
    int k=subGroup.size();
    for(int i=0; i<k; i++){
        float maxQ=0; int idx = 0;
        for(int j=0; j<subGroup.size(); j++)
            if(subGroup[j].q>maxQ){
                maxQ = subGroup[j].q;
                idx = j;
            }
        subGroupSorted.push_back(subGroup[idx]);
        subGroup.erase(subGroup.begin()+idx);
    }

    ofstream myfile("OUTPUT");
    // print header
    myfile <<setw(7)<<left<<"k"<<setw(6)<<"p, p0="<<setw(9)<<((int)(p0*100))/100.f<<setw(5)<<"q, F="<<setw(10)<<db.F<<setw(16)<<"V, V_a=2.576"<<"hypothesis, n="<<db.n<<endl;
    // sub-group loop
    for(int i=0; i<subGroupSorted.size(); i++){
        float V = (subGroupSorted[i].p-p0) / sqrt(p0*(1-p0)/subGroupSorted[i].g/db.m); // z-transform
        // print properties
        myfile <<setw(7)<<left<<i+1<<setw(15)<<subGroupSorted[i].p<<setw(15)<<subGroupSorted[i].q<<setw(2)<<((max(V,-V)>2.576) ? "  " : "X ")<<setw(14)<<V;
        // print hypothesis
        for(int j=0; j<subGroupSorted[i].groupHypo.size(); j++)
            myfile << setw(2) << right << subGroupSorted[i].groupHypo[j].attribute << "=" << subGroupSorted[i].groupHypo[j].value << " ";
        myfile << endl;
    }
    myfile.close();
    cout<<"Sub-group saved in ./OUTPUT"<<endl<<endl;
}
