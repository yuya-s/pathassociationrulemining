
#ifndef PATHASSOCIATION_UTILITY_H
#define PATHASSOCIATION_UTILITY_H

#include "graph.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <random>
#include <algorithm>
#include <map>
#include <set>
#include <list>
#include <omp.h>
#include <mutex>
#include <google/dense_hash_map>


using namespace std;

class Result {
public:

    string outputFile;
    string inputFile;

    double minsup;
    double miningTime;
    int maxlength;
    double totalminingTime;
    double attributeTime;
    double lengthoneTime;
    double lengthLTime;
    double kleeneTime;
    double combineTime;
    int pathNum;
    int patternNum;
    int coreNum;
    int efficientFlag;
    int extentionFlag;
    int stratifiedsampleFlag;
    int partitionmode;
    double approximatecandidaterate;
    double approximatesamplingrate;

//    int querynum;
//    double querySumTime;
//    double maintenanceAddEdgeTime;
//    double maintenanceDeleteEdgeTime;
//    double maintenanceAddWorkloadTime;
//    double maintenanceDeleteWorkloadTime;
//    int maintenanceEdgeNum;
//    int maintenanceWorkloadNum;

    void clear(){
        minsup=0;
        miningTime=0;
        totalminingTime=0;
        pathNum=0;
        efficientFlag=0;
    }
    void Output() {
        std::ofstream fout(outputFile, ios::app);
        fout << inputFile<<"\t" << efficientFlag <<"\t"<<  extentionFlag <<"\t"<<  stratifiedsampleFlag <<"\t"<<  approximatecandidaterate <<"\t"<<  approximatesamplingrate <<"\t"<< coreNum<< "\t" << minsup<<"\t"<<maxlength<< "\t"<<partitionmode<< "\t" << attributeTime <<"\t"<<lengthoneTime <<"\t"<<lengthLTime <<"\t"<<kleeneTime <<"\t"<<combineTime   << "\t" << miningTime<< "\t" << pathNum<< "\t" << patternNum<<  std::endl;
    }


//    void OutputMaintenance(){
//        std::ofstream fout(outputFile+"_maintenance", ios::app);
//        //fout << inputFile << "\t" << k  << "\t"<<maintenanceEdgeNum<<"\t"<<maintenanceWorkloadNum<<"\t" << maintenanceAddEdgeTime/(double)maintenanceEdgeNum << "\t" << maintenanceDeleteEdgeTime/(double)maintenanceEdgeNum<< "\t" << maintenanceAddWorkloadTime/(double)maintenanceWorkloadNum << "\t" << maintenanceDeleteWorkloadTime/(double)maintenanceWorkloadNum<<std::endl;
//    }

};

inline std::vector<std::string> split(std::string str, char del) {
    int first = 0;
    int last = str.find_first_of(del);

    std::vector<std::string> result;

    while (first < str.size()) {
        std::string subStr(str, first, last - first);

        result.push_back(subStr);

        first = last + 1;
        last = str.find_first_of(del, first);

        if (last == std::string::npos) {
            last = str.size();
        }
    }

    return result;
}


static bool Contains(const vector<int>& vec, int element)
{
	return std::find(vec.begin(), vec.end(), element) != vec.end();
}

static void ShowVector(const vector<int>& vec)
{
    for(int i=0; i< vec.size();i++){
        cout<<vec[i]<<" ";
    }
    cout<<endl;
}



static int Contains_ReturnIndex(const vector<int>& vec, int element)
{
    for(unsigned int i=0; i < vec.size();i++){
        if(vec[i]==element)return i;
    }
    return -1;
}

static bool ContainsVector(const vector<vector<int>>& vecset, vector<int> vec)
{
	return std::find(vecset.begin(), vecset.end(), vec) != vecset.end();
}

static bool Overlap(const vector<int>& vec1, const vector<int>& vec2)
{
    for(unsigned int i=0; i < vec1.size();i++) {
        for (unsigned int j = 0; j < vec2.size(); j++) {
            if (vec1[i] == vec2[j])return true;
        }
    }
    return false;
}

template<class T> bool vectorExist (vector<T> c, T item)
{
    return (std::find(c.begin(), c.end(), item) != c.end());
}

template<class T> vector<T> vectorUnion (vector<T> a, vector<T> b)
{
    vector<T> c;

    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());

    auto i = a.begin();
    auto j = b.begin();

    while (i != a.end() || j != b.end())
    {
        if (j == b.end() || *i < *j)
        {
            if(!vectorExist(c,*i)) c.push_back(*i);
            i++;
        }
        else
        {
            if(!vectorExist(c,*j)) c.push_back(*j);
            j++;
        }
    }

    return c;
}

static bool AllOverlap(const vector<int>& vec1, const vector<int>& vec2)
{
    if(vec1.size()>vec2.size()){

        for(unsigned int i=0; i < vec2.size();i++) {
            if(!Contains(vec1,vec2[i]))return false;
        }

    }
    else{
        for(unsigned int i=0; i < vec1.size();i++) {
            if(!Contains(vec2,vec1[i]))return false;
        }
    }

    return true;
}

static vector<int> intersection(vector<int> &v1, vector<int> &v2){
    
    std::vector<int> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),back_inserter(v3));
    return v3;
}


static bool Vec1AllOverlap(const vector<int>& vec1, const vector<int>& vec2)
{
    if(vec1.size()>=vec2.size()){

        for(unsigned int i=0; i < vec2.size();i++) {
            if(!Contains(vec1,vec2[i]))return false;
        }

    }
    else{
        return false;
    }

    return true;
}


static bool Merge(const vector<int>& vec1, const vector<int>& vec2)
{
    for(unsigned int i=0; i < vec1.size();i++) {
        for (unsigned int j = 0; j < vec2.size(); j++) {
            if (vec1[i] == vec2[j])return true;
        }
    }
    return false;
}

static void Remove(vector<int>* vec, int element)
{
	vec->erase(std::remove(vec->begin(), vec->end(), element), vec->end());
}

template<class T>
void combination(const vector<T>& seed, int target_size, vector<vector<T>>& combinations) {
    vector<int> indices(target_size);
    const int seed_size = seed.size();
    int start_index = 0;
    int size = 0;

    while (size >= 0) {
        for (int i = start_index; i < seed_size; ++i) {
            indices[size++] = i;
            if (size == target_size) {
                vector<T> comb(target_size);
                for (int x = 0; x < target_size; ++x) {
                    comb[x] = seed[indices[x]];
                }
                combinations.push_back(comb);
                break;
            }
        }
        --size;
        if (size < 0) break;
        start_index = indices[size] + 1;
    }
}

#endif //PATHASSOCIATION_UTILITY_H
