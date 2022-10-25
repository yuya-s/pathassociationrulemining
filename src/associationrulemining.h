
#ifndef PATHASSOCIATION_ASSOCIATIONRULEMINING_H
#define PATHASSOCIATION_ASSOCIATIONRULEMINING_H

#include "graph.h"

#define RANDOM 0
#define AVERAGESCORE 1
#define AVERAGESCOREVERTEXNUMBALANCE 2


class Path{

public:
    int pathidentifier;
    int count;
    vector<vector<int>> vertextuples;
    vector<int> edgelabels;
    bool kleenestar;
    vector<int> verticalExtendingPathidentifiers;
    vector<int> horizontalExtendingPathidentifiers;

    void clear(){

        pathidentifier=-1;
        count=0;
        vertextuples.clear();
        edgelabels.clear();
        kleenestar=false;
    }

    bool operator==(const Path& p){
        if(kleenestar!=p.kleenestar)return false;
        if(edgelabels.size()!=p.edgelabels.size())return false;
        for(int i = 0; i < edgelabels.size();i++){
            if(edgelabels[i]!=p.edgelabels[i])return false;
        }
        for(int i = 0; i < vertextuples.size();i++){
            if(vertextuples[i]!=p.vertextuples[i])return false;
        }
        return true;
//        return(vertextuples==p.vertextuples&&edgelabels==p.edgelabels&&kleenestar==p.kleenestar);
    }
    bool operator>(const Path& p){
        if(!kleenestar&&p.kleenestar)return true;
        else if (kleenestar&&!p.kleenestar)return false;

        if(edgelabels.size()<p.edgelabels.size())return true;
        else if(edgelabels.size()>p.edgelabels.size())return false;

        for(int i = 0; i < edgelabels.size();i++){
            if(edgelabels[i]<p.edgelabels[i])return true;
            else if(edgelabels[i]>p.edgelabels[i])return false;
        }
        for(int i = 0; i < vertextuples.size();i++){
            if(vertextuples[i]<p.vertextuples[i])return true;
            else if(vertextuples[i]>p.vertextuples[i])return false;
        }
        return false;
//        return(vertextuples==p.vertextuples&&edgelabels==p.edgelabels&&kleenestar==p.kleenestar);
    }


    void Show(){
         cout<<pathidentifier<<": ";
            for(int for_vertex=0;for_vertex<vertextuples.size();for_vertex++){

                cout<<"( ";
                for(int for_tuple=0;for_tuple<vertextuples[for_vertex].size();for_tuple++){

                    cout<<vertextuples[for_vertex][for_tuple]<<" ";

                }
                cout<<") ";

                if(for_vertex!=vertextuples.size()-1){
                    if(kleenestar)cout<<edgelabels[for_vertex]<<"* ";
                    else cout<<edgelabels[for_vertex]<<" ";
                }
            }
            cout<<"count="<<count<<endl;
    }
};

class PatternPath{

public:
    int source;
    int target;
    int pathidentifier;
};

class Pattern{

public:
    int patternidentifier;
    int count;
    vector<int> pathidentifiers;
    bool duplicateFlag;

    vector<vector<int>> vertextuples;
    vector<PatternPath> patternedges;

    double confidence1;
    double confidence2;

    double lift;

    Pattern()
	{
        patternidentifier = 0;
		count = 0;
		pathidentifiers.clear();
        duplicateFlag = false;
	}

    bool operator==(const Pattern& p){
        if(pathidentifiers[0]==p.pathidentifiers[0]&&pathidentifiers[1]==p.pathidentifiers[1])return true;
        if(pathidentifiers[1]==p.pathidentifiers[0]&&pathidentifiers[0]==p.pathidentifiers[1])return true;
        return false;
//        return(vertextuples==p.vertextuples&&edgelabels==p.edgelabels&&kleenestar==p.kleenestar);
    }



};

class MatchedVertex{

public:
    int vertexid;

    vector<vector<int>> pathidentifiers;
    vector<int> kleenepaths;
    vector<vector<vector<int>>> targets;
    vector<int> patternidentifiers;
    vector<vector<int>> kleenetagets;
    vector<int> targetidentifiers;

    int number_of_frequentedges;
    int number_of_frequentattributes;

};

class AssociationRules{

public:
    unsigned int maximumpathlength;
    unsigned int minsupport;
    unsigned int corenum;

    vector<vector<Path>> paths;
    vector<Pattern> patterns;
    vector<MatchedVertex> matchedvertices;
    int globalpathidcount;
    int globaltargetcount;

    vector<vector<int>> parallelvertexidsets;

    vector<vector<bool>> matchedattributesets_vertex;
    vector<vector<int>> matchedvertexs_attributeset;
    vector<vector<int>>candidateedgelabels;
    vector<int> frequentattributes;
    vector<vector<vector<pair<int,vector<int>>>>> targettuples_edges;

    google::dense_hash_map<int, pair<int,int>> pathidentifier2pathsindex;

    int maximumsizemultipleattributes;

//    double minconfidence;
    bool s_join;
    bool s_t_join;
    bool kleenestar;

    struct timespec starttime;
    double timelimit;

    void CheckTime(){
        struct timespec currenttime;
        clock_gettime(CLOCK_REALTIME, &currenttime);
        double elapsedtime= (currenttime.tv_sec - starttime.tv_sec) + (currenttime.tv_nsec - starttime.tv_nsec) * pow(10, -9);
        //cout<<elapsedtime<<"/"<<timelimit*3600<<endl;
        if(timelimit*3600<elapsedtime){
            cout<< "over " <<timelimit<< "hours"<<endl;
            exit(1);
        }
    };

    void FindPattern(OriginalGraph);

    bool VertexAttributeMatching(vector<int> &graphvertexattributes_, vector<int> &rulevertexattributes_);
    bool CheckPathDominance(Path& path1, Path& path2);//if path1 (resp. path2) dominates path2 (resp. path1), return true;
    bool CheckVertexTupleDominance(vector<int>& vertex1, vector<int>& vertex2); //if vertex1 dominates vertex2, return true;


//    void FindFrequentSingleAttribute(OriginalGraph&);
//    void FindFrequentMultipleAttribute(OriginalGraph&);
//    void FindFrequentLengthOnePath(OriginalGraph&, vector<Path>&, int);
//    void FindFrequentLengthLPath(OriginalGraph&, int length, int maxattributes);
//    void FindFrequentKleenePath(OriginalGraph&);
//    void CombinePaths(OriginalGraph);

    //void FindCandidateEdgeLabel(OriginalGraph&);
    void ShowAllPaths(string);
    void ShowAllPatterns(string);
    void SimpleShowAllPatterns(string);

    Path MergeTwoPathSourceAndTarget(Path&, Path&, int);
    Path MergeTwoOneLengthPathSource(Path&, Path&, int);
    Path MergeTwoOneLengthPathTarget(Path&, Path&, int);
    Path MergeTwoPath(Path&, Path&, int);
    Path MergeTwoPathGivenTuples(Path& path1, Path& path2, int numberOfCandidateAttribute, vector<int>& mergetuple);


    Path MergeTwoLLengthPathSource(Path&, Path&, int length, int numatt);
    Path MergeTwoLLengthPathTarget(Path&, Path&, int length, int numatt);

//    void EfficientFindFrequentMultipleAttribute(OriginalGraph&);
//    void EfficientFindFrequentLengthOnePath(OriginalGraph &);
//    void EfficientFindFrequentLengthLPath(OriginalGraph &, int);
//    void EfficientFindFrequentKleenePath(OriginalGraph&);
//    void EfficientCombinePaths(OriginalGraph);




    void ParallelFindFrequentSingleAttribute(OriginalGraph&);
    void ParallelFindFrequentMultipleAttribute(OriginalGraph&);
    void ParallelFindFrequentLengthOnePath(OriginalGraph &, vector<Path>&, int);
    void ParallelFindFrequentLengthLPath(OriginalGraph&, int length, vector<Path>&, int maxattributes);
    void ParallelFindFrequentKleenePath(OriginalGraph&);

    void ParallelEfficientFindFrequentMultipleAttribute(OriginalGraph &);
    void ParallelEfficientFindFrequentLengthOnePath(OriginalGraph &);
    void ParallelEfficientFindFrequentLengthLPath(OriginalGraph &, int);
    void ParallelEfficientFindFrequentKleenePath(OriginalGraph&);
    bool CandidatePathMatching(Path& path, vector<int>& pathidentifiers, vector<int>& mergetuples, int position, int start, int end, int length);

    void ParallelEfficientCombinePaths(OriginalGraph);
    void ParallelCombinePaths(OriginalGraph);

    void ParallelVertexAttributeMatching(OriginalGraph& graph, vector<int>& matchedvertexids, vector<int>& vertexidset, int attribute, vector<int> & countinverseedgelabels);
    vector<int> ParallelVertexMultipleAttributeMatching(OriginalGraph& graph, vector<int>& vertexidset, vector<int>& pathidset);
    vector<int> ParallelVertexMatching(OriginalGraph &, vector<pair<int,int>>& vertexidset, vector<vector<int>>& taregetslist, vector<vector<int>>& nexttarget, int path0id);
    vector<int> ParallelVertexMatching_forEfficient(OriginalGraph &, vector<int>& vertexidset, vector<int>& pathidset, int length, vector<vector<int>>& nexttarget);
    vector<int> PatternVertexMatching_forEfficient(OriginalGraph &, vector<int>& vertexidset, int length1, int length2, int pathid1, int pathid2);
    void ParallelKleenePathMatching_perEdgeLabel(OriginalGraph& graph, int edgelabel);

    vector<vector<int>> ParallelFindTarget(OriginalGraph &, vector<pair<int,int>>& vertexidset, int edgelabel, int length);
    vector<pair<int,int>> ParallelPathContainment(vector<int>&, int length, int pathidentifier);

    void VerticesPartitioning(OriginalGraph & graph, int mode, bool initial);

    void MetricsComputation(OriginalGraph& graph);
    void DuplicateRemove();

    //void FindFrequentLengthLPath_test(OriginalGraph&, int);
};

#endif //PATHASSOCIATION_ASSOCIATIONRULEMINING_H
