

#include "associationrulemining.h"


extern bool EfficientFlag;
extern bool ExtentionFlag;
extern int PartitionMode;
extern double ApproximateCandidateRate;
extern double ApproximateSamplingRate;
extern bool StratifiedSmaplingFlag;
extern bool EmptysourceFlag;
extern bool NonDominanceOutputFlag;
extern int minimumRulePathLength;
extern double attributetime;
extern double lengthonetime;
extern double lengthLtime;
extern double kleenetime;
extern double combinetime;


void AssociationRules::FindPattern(OriginalGraph graph){

    struct timespec frequentattribute_startTime, frequentattribute_endTime, frequentpath__startTime, frequentpath_endTime, frequentkleenepath_startTime, frequentkleenepath_endTime;
    struct timespec frequentpattern_startTime, frequentpattern_endTime;
    struct timespec temp1_startTime, temp1_endTime;
    struct timespec temp2_startTime, temp2_endTime;
    struct timespec temp3_startTime, temp3_endTime;
    struct timespec temp4_startTime, temp4_endTime;
    globalpathidcount=0;
    globaltargetcount=0;
    clock_gettime(CLOCK_REALTIME, &starttime);
    //initialize//
    paths.resize(maximumpathlength+1);
    for(int i=0;i<maximumpathlength+1;i++)paths[i].clear();
    patterns.clear();
    matchedvertices.resize(graph.number_of_vertices);
    matchedattributesets_vertex.resize(graph.number_of_vertices);
    frequentattributes.clear();

    for(int i = 0; i<graph.number_of_vertices;i++){

        matchedvertices[i].vertexid=i;

        matchedvertices[i].pathidentifiers.resize(maximumpathlength+1);
        matchedvertices[i].targets.resize(maximumpathlength+1);
        for(int j=0;j<maximumpathlength+1;j++){
            matchedvertices[i].pathidentifiers[j].clear();
            matchedvertices[i].targets[j].clear();
        }

        matchedvertices[i].patternidentifiers.clear();
    }

    candidateedgelabels.resize(maximumpathlength);
    targettuples_edges.resize(maximumpathlength);
    for(int i=0;i<maximumpathlength;i++) {
        candidateedgelabels.clear();
        targettuples_edges[i].resize(graph.cardinality_of_edgelabels);
        for (int j = 0; j < graph.cardinality_of_edgelabels; j++)targettuples_edges[i][j].clear();
    }

    pathidentifier2pathsindex.set_empty_key(-1);

    //cout<<"efficient flag "<<EfficientFlag<<endl;

    VerticesPartitioning(graph, RANDOM, true);


    //**Length 0 Frequent Paths i.e., Vertex attributes **//
    {
        double time_single,time_multiple;
        int size_single, size_multiple;
        cout << "START finding frequent length 0 frequent paths" << endl;
        clock_gettime(CLOCK_REALTIME, &frequentattribute_startTime);
        //each attribute//
        cout << "---single attrebutes" << endl;
        ParallelFindFrequentSingleAttribute(graph);

        clock_gettime(CLOCK_REALTIME, &frequentattribute_endTime);
        size_single=paths[0].size();
        time_single = (frequentattribute_endTime.tv_sec - frequentattribute_startTime.tv_sec) + (frequentattribute_endTime.tv_nsec - frequentattribute_startTime.tv_nsec) * pow(10, -9);

        vector<int> targetnum(maximumpathlength,0);
        for(int i=0;i<maximumpathlength;i++){
            for(int j=0;j<graph.cardinality_of_edgelabels;j++){
                targetnum[i]+= targettuples_edges[i][j].size();
            }
        }
        cout << "---DONE: single attribute size= "<<size_single<<", target size = (";
        for(int i=0;i<maximumpathlength;i++)cout<<" "<<targetnum[i]<<",";
        cout<<") edge label size = (";
        for(int i=0;i<maximumpathlength;i++)cout<<" "<< candidateedgelabels[i].size()<<",";
        cout <<"), time=" << time_single << endl;

        if(paths[0].size()==0)return; // No paths mean no need to search further

        clock_gettime(CLOCK_REALTIME, &frequentattribute_startTime);
        cout << "---multiple attributes" << endl;

        if(EfficientFlag)ParallelEfficientFindFrequentMultipleAttribute(graph);
        else ParallelFindFrequentMultipleAttribute(graph);

        clock_gettime(CLOCK_REALTIME, &frequentattribute_endTime);
        time_multiple = (frequentattribute_endTime.tv_sec - frequentattribute_startTime.tv_sec) + (frequentattribute_endTime.tv_nsec - frequentattribute_startTime.tv_nsec) * pow(10, -9);
        size_multiple=paths[0].size() - size_single;
        cout << "---DONE: size= "<<size_multiple<<", time=" << time_multiple << endl;
        cout << "DONE: path0 size=" << paths[0].size() << ", time=" << time_single + time_multiple << endl;

        attributetime=time_single + time_multiple;
    }

    if(PartitionMode!=RANDOM || ApproximateSamplingRate!=1.0){
        cout<<"Reparitioning :" << PartitionMode<<endl;
        VerticesPartitioning(graph, PartitionMode, false);
    }

    for(int i=0;i<graph.number_of_vertices;i++){
        matchedattributesets_vertex[i].resize(paths[0].size());
        for(int j=0;j<paths[0].size();j++)matchedattributesets_vertex[i][j]=false;
    }
    for(int for_path=0;for_path<paths[0].size();for_path++){
        for(auto& vertex: matchedvertexs_attributeset[for_path]){
            matchedattributesets_vertex[vertex][for_path]=true;
        }
    }

    clock_gettime(CLOCK_REALTIME, &frequentpath__startTime);

    {
        cout << "START finding length 1 frequent paths" << endl;

        if(EfficientFlag)ParallelEfficientFindFrequentLengthOnePath(graph);
        else ParallelFindFrequentLengthOnePath(graph,paths[0],graph.number_of_vertices);

        cout << "DONE: path 1 size= " << paths[1].size() << ", time= ";
        clock_gettime(CLOCK_REALTIME, &frequentpath_endTime);
        lengthonetime=(frequentpath_endTime.tv_sec - frequentpath__startTime.tv_sec) + (frequentpath_endTime.tv_nsec - frequentpath__startTime.tv_nsec) * pow(10, -9);
        cout << lengthonetime<<endl;

        cout << "START finding length L frequent paths" << endl;
        for (int for_length = 2; for_length < maximumpathlength + 1; for_length++) {
            if(EfficientFlag)ParallelEfficientFindFrequentLengthLPath(graph, for_length);
            else ParallelFindFrequentLengthLPath(graph, for_length,paths[for_length-1],graph.number_of_vertices);
            //else ParallelFindFrequentLengthLPath(graph, for_length,paths[for_length-1],graph.number_of_vertices);
            //else FindFrequentLengthLPath(graph, for_length,graph.number_of_vertices);

            clock_gettime(CLOCK_REALTIME, &frequentpath_endTime);
            cout << "DONE: path" << for_length << ", size= " << paths[for_length].size() << ", time= ";
            lengthLtime=(frequentpath_endTime.tv_sec - frequentpath__startTime.tv_sec) + (frequentpath_endTime.tv_nsec - frequentpath__startTime.tv_nsec) * pow(10, -9);
            cout << lengthLtime << endl;
        }
    }

    cout<<"START: kleene path fiding"<<endl;
    clock_gettime(CLOCK_REALTIME, &frequentkleenepath_startTime);

    if(maximumpathlength>1 && kleenestar)
    {
        if (EfficientFlag)ParallelEfficientFindFrequentKleenePath(graph);
        //else if (EfficientFlag)EfficientFindFrequentKleenePath(graph);
        else ParallelFindFrequentKleenePath(graph);
        //else FindFrequentKleenePath(graph);
    }
    clock_gettime(CLOCK_REALTIME, &frequentkleenepath_endTime);
    cout<<"DONE: "<< paths[1].size()<< ", time= ";
    kleenetime=(frequentkleenepath_endTime.tv_sec - frequentkleenepath_startTime.tv_sec) + (frequentkleenepath_endTime.tv_nsec - frequentkleenepath_startTime.tv_nsec) * pow(10, -9);

    cout<<kleenetime<<endl;

    cout<<"START: finding frequent patterns"<<endl;
    clock_gettime(CLOCK_REALTIME, &frequentpattern_startTime);

    if(EfficientFlag)ParallelEfficientCombinePaths(graph);
    //else if(EfficientFlag)EfficientCombinePaths(graph);
    else ParallelCombinePaths(graph);
    //else CombinePaths(graph);

    clock_gettime(CLOCK_REALTIME, &frequentpattern_endTime);
    cout<<"DONE: time= ";
    combinetime=(frequentpattern_endTime.tv_sec - frequentpattern_startTime.tv_sec) + (frequentpattern_endTime.tv_nsec - frequentpattern_startTime.tv_nsec) * pow(10, -9);
    cout<<combinetime<<endl;


//    cout<<"finish"<<endl;

}


void AssociationRules::ParallelFindFrequentMultipleAttribute(OriginalGraph& graph){

    struct timespec temp1_startTime, temp1_endTime;
    double temp1time;
    clock_gettime(CLOCK_REALTIME, &temp1_startTime);

    int start=0;
    int end=paths[0].size();
    int pathlength=1;
    int maximumsizecombination=paths[0].size();
    int numberOfCandidateAttribute=1;
    vector<vector<int>>tempcandidatetuplesvector;
    vector<vector<int>>candidatetuplesvector;

    while(end!=0) {

        tempcandidatetuplesvector.clear();
        candidatetuplesvector.clear();
        //Candidate generation
        for (int for_path1 = start; for_path1 < end - 1; for_path1++) {
            if(EmptysourceFlag&&paths[0][for_path1].vertextuples[0].empty())continue;

            for (int for_path2 = for_path1 + 1; for_path2 < end; for_path2++) {

                int numberOfMathcedAttribute=0;
                bool diffFlag=true;
                for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
                    diffFlag=false;
                    if(paths[0][for_path1].vertextuples[0][for_att] != paths[0][for_path2].vertextuples[0][for_att]) {
                        numberOfCandidateAttribute=for_att;
                        diffFlag=true;
                        break;
                    }
                }

                if(numberOfMathcedAttribute == numberOfCandidateAttribute-1&&diffFlag){
                    vector<int> candidatetuple;
                    for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
                        candidatetuple.push_back(paths[0][for_path1].vertextuples[0][for_att]);
                    }

                    int att1 = paths[0][for_path1].vertextuples[0][numberOfCandidateAttribute-1];
                    int att2 = paths[0][for_path2].vertextuples[0][numberOfCandidateAttribute-1];
                    if(att1>att2)swap(att1,att2);

                    candidatetuple.push_back(att1),candidatetuple.push_back(att2);
                    tempcandidatetuplesvector.push_back(candidatetuple);
                }
            }
        }

        for(auto& tuple: tempcandidatetuplesvector){

            bool containFlag=false;
            bool insertFlag=true;
            for(int for_att=0;for_att<tuple.size();for_att++){

                vector<int> tmptuple = tuple;
                tmptuple.erase(tmptuple.begin()+for_att);

                for (int for_path = start; for_path < end; for_path++) {

                    if(paths[0][for_path].vertextuples[0]==tmptuple){
                        containFlag=true;
                        break;
                    }
                }

                if(!containFlag){
                    insertFlag=false;
                    break;
                }
            }
            if(insertFlag){
                candidatetuplesvector.push_back(tuple);
            }
        }

        if(candidatetuplesvector.size()==0)break;

        //Find frequent Set//
        for(auto& rulevertextuple: candidatetuplesvector){
            int count = 0;
//            vector<int> matchedvertexids;
//            matchedvertexids.clear();
            vector<vector<int>> matchedvertexids;
            matchedvertexids.resize(corenum);


            #pragma omp parallel for
            for(int for_core=0;for_core<corenum;for_core++) {
                 matchedvertexids[for_core].clear();
                 matchedvertexids[for_core]=ParallelVertexMultipleAttributeMatching(graph, parallelvertexidsets[for_core],  rulevertextuple);
            }

            for(int for_core=0;for_core<corenum;for_core++) {
                count+=matchedvertexids[for_core].size();
            }

//            for (auto vertex: graph.vertex_list) {
//
//                bool matchFlag;
//
//                matchFlag = VertexAttributeMatching(vertex.vertextuple, rulevertextuple);
//
//                if(matchFlag){
//                    count++;
//                    matchedvertexids.push_back(vertex.vertexid);
//                }
//            }

            if(count>minsupport){

                Path temppath;
                temppath.pathidentifier=globalpathidcount;
                pathidentifier2pathsindex.insert({globalpathidcount,make_pair(0,paths[0].size())});
                globalpathidcount++;
                temppath.count=count;
                temppath.vertextuples.push_back(rulevertextuple);
                temppath.kleenestar=false;

                paths[0].push_back(temppath);

                #pragma omp parallel for
                for(int for_core=0;for_core<corenum;for_core++) {
                    for (auto matchedvertexid: matchedvertexids[for_core]) {

                        matchedvertices[matchedvertexid].pathidentifiers[0].push_back(temppath.pathidentifier);
                        vector<int> temptarget;
                        temptarget.push_back(matchedvertexid);
                        matchedvertices[matchedvertexid].targets[0].push_back(temptarget);
                    }
                }

//                for(auto matchedvertexid: matchedvertexids){
//
//                    matchedvertices[matchedvertexid].pathidentifiers[0].push_back(temppath.pathidentifier);
//                    vector<int> temptarget;
//                    temptarget.push_back(matchedvertexid);
//                    matchedvertices[matchedvertexid].targets[0].push_back(temptarget);
//                    //matchedvertices[matchedvertexid].patternidentifiers.push_back(temppattern.patternidentifier);
//                }
                vector<int> allmatchedvertexids;
                allmatchedvertexids.clear();
                for(int for_core=0;for_core<corenum;for_core++) {
                    for(auto& vertexid: matchedvertexids[for_core]) {
                        allmatchedvertexids.push_back(vertexid);
                    }
                }
                matchedvertexs_attributeset.push_back(allmatchedvertexids);
            }
        }
        start=end;
        end=paths[0].size();
        numberOfCandidateAttribute++;
        if(start==end)break;
    }

    clock_gettime(CLOCK_REALTIME, &temp1_endTime);
    temp1time=(temp1_endTime.tv_sec - temp1_startTime.tv_sec) + (temp1_endTime.tv_nsec - temp1_startTime.tv_nsec) * pow(10, -9);
    cout << "multi attribute temp1 = "<<temp1time<<endl;

    for(int for_edgelabel=0;for_edgelabel<graph.cardinality_of_edgelabels;for_edgelabel++) {

        int start=0;
        int end = targettuples_edges[maximumpathlength-1][for_edgelabel].size();
        numberOfCandidateAttribute=1;

        while (end != 0) {

            tempcandidatetuplesvector.clear();
            candidatetuplesvector.clear();
            //Candidate generation
            for (int for_tuple1 = start; for_tuple1 < end - 1; for_tuple1++) {

                for (int for_tuple2 = for_tuple1 + 1; for_tuple2 < end; for_tuple2++) {

                    int numberOfMathcedAttribute = 0;
                    bool diffFlag = true;
                    for (int for_att = 0; for_att < numberOfCandidateAttribute - 1; for_att++) {
                        diffFlag = false;
                        if (targettuples_edges[maximumpathlength-1][for_edgelabel][for_tuple1].second[for_att] != targettuples_edges[maximumpathlength-1][for_edgelabel][for_tuple2].second[for_att]) {
                            numberOfCandidateAttribute = for_att;
                            diffFlag = true;
                            break;
                        }
                    }

                    if (numberOfMathcedAttribute == numberOfCandidateAttribute - 1 && diffFlag) {
                        vector<int> candidatetuple;
                        for (int for_att = 0; for_att < numberOfCandidateAttribute - 1; for_att++) {
                            candidatetuple.push_back(targettuples_edges[maximumpathlength-1][for_edgelabel][for_tuple1].second[for_att]);
                        }

                        int att1 = targettuples_edges[maximumpathlength-1][for_edgelabel][for_tuple1].second[numberOfCandidateAttribute - 1];
                        int att2 = targettuples_edges[maximumpathlength-1][for_edgelabel][for_tuple2].second[numberOfCandidateAttribute - 1];
                        if (att1 > att2)swap(att1, att2);

                        candidatetuple.push_back(att1), candidatetuple.push_back(att2);
                        tempcandidatetuplesvector.push_back(candidatetuple);
                    }
                }
            }
            //cout<<tempcandidatetuplesvector.size()<<endl;

            for (auto &tuple: tempcandidatetuplesvector) {

                bool containFlag = false;
                bool insertFlag = true;
                for (int for_att = 0; for_att < tuple.size(); for_att++) {

                    vector<int> tmptuple = tuple;
                    tmptuple.erase(tmptuple.begin() + for_att);

                    for (int for_tuple = start; for_tuple < end; for_tuple++) {

                        if (targettuples_edges[maximumpathlength-1][for_edgelabel][for_tuple].second == tmptuple) {
                            containFlag = true;
                            break;
                        }
                    }

                    if (!containFlag) {
                        insertFlag = false;
                        break;
                    }
                }
                if (insertFlag) {
                    candidatetuplesvector.push_back(tuple);
                }
            }

            if (candidatetuplesvector.size() == 0)break;

            //Find frequent Set//
            //cout<<"size="<<candidatetuplesvector.size()<<endl;
            for (auto &rulevertextuple: candidatetuplesvector) {
                int count = 0;
                vector<int> matchedvertexids;
                matchedvertexids.clear();
                for (auto vertex: graph.vertex_list) {

                    bool matchFlag;

                    matchFlag = VertexAttributeMatching(vertex.vertextuple, rulevertextuple);

                    if (matchFlag) {
                        count+=graph.vertex_inverseconnectingcount_edgelabels[vertex.vertexid][for_edgelabel];
                        matchedvertexids.push_back(vertex.vertexid);
                    }
                }
                for(int for_length=0;for_length<maximumpathlength;for_length++){

                    if(count!=0){

                        for(int for_length2=for_length;for_length2<maximumpathlength;for_length2++){
                            targettuples_edges[for_length2][for_edgelabel].push_back(make_pair(globaltargetcount,rulevertextuple));
                        }
                        for(auto& matchedvertexid: matchedvertexids){
                            matchedvertices[matchedvertexid].targetidentifiers.push_back(globaltargetcount);
                        }
                        globaltargetcount++;
                        break;
                    }

                }
            }
            start = end;
            end = targettuples_edges[maximumpathlength-1][for_edgelabel].size();
            numberOfCandidateAttribute++;
            if (start == end)break;
        }
    }

    maximumsizemultipleattributes=paths[0].back().vertextuples[0].size();
}


void AssociationRules::ParallelEfficientFindFrequentMultipleAttribute(OriginalGraph& graph){


    int start=0;
    int end=paths[0].size();
    int pathlength=1;
    int maximumsizecombination=paths[0].size();
    int numberOfCandidateAttribute=1;
    vector<pair<vector<int>,vector<int>>>tempcandidatetuplesvector;
    vector<pair<vector<int>,vector<int>>>candidatetuplesvector;

    while(end!=0) {

        tempcandidatetuplesvector.clear();
        candidatetuplesvector.clear();
        //Candidate generation
        for (int for_path1 = start; for_path1 < end - 1; for_path1++) {
            if(EmptysourceFlag&&paths[0][for_path1].vertextuples[0].empty())continue;
            for (int for_path2 = for_path1 + 1; for_path2 < end; for_path2++) {

                int numberOfMathcedAttribute=0;
                bool diffFlag=true;
                for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
                    diffFlag=false;
                    if(paths[0][for_path1].vertextuples[0][for_att] != paths[0][for_path2].vertextuples[0][for_att]) {
                        numberOfCandidateAttribute=for_att;
                        diffFlag=true;
                        break;
                    }
                }

                if(numberOfMathcedAttribute == numberOfCandidateAttribute-1&&diffFlag){
                    vector<int> candidatetuple;
                    vector<int> candidatepathidentifiers;
                    for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
                        candidatetuple.push_back(paths[0][for_path1].vertextuples[0][for_att]);
                    }

                    int att1 = paths[0][for_path1].vertextuples[0][numberOfCandidateAttribute-1];
                    int att2 = paths[0][for_path2].vertextuples[0][numberOfCandidateAttribute-1];
                    if(att1>att2)swap(att1,att2);

                    candidatetuple.push_back(att1),candidatetuple.push_back(att2);
                    candidatepathidentifiers.push_back(paths[0][for_path1].pathidentifier);
                    candidatepathidentifiers.push_back(paths[0][for_path2].pathidentifier);
                    tempcandidatetuplesvector.push_back(make_pair(candidatetuple,candidatepathidentifiers));
                }
            }
        }

        for(auto& tuple: tempcandidatetuplesvector){

            bool containFlag=false;
            bool insertFlag=true;
            for(int for_att=0;for_att<tuple.first.size();for_att++){

                vector<int> tmptuple = tuple.first;
                tmptuple.erase(tmptuple.begin()+for_att);

                for (int for_path = start; for_path < end; for_path++) {

                    if(paths[0][for_path].vertextuples[0]==tmptuple){
                        containFlag=true;
                        break;
                    }
                }

                if(!containFlag){
                    insertFlag=false;
                    break;
                }
            }
            if(insertFlag){
                candidatetuplesvector.push_back(tuple);
            }
        }

        if(candidatetuplesvector.size()==0)break;

        //Find frequent Set//
        for(auto& rulevertextuple: candidatetuplesvector){
            int count = 0;

            vector<vector<int>> matchedvertexids;
            matchedvertexids.resize(corenum);


            #pragma omp parallel for
            for(int for_core=0;for_core<corenum;for_core++) {
                 matchedvertexids[for_core].clear();
                 matchedvertexids[for_core]=ParallelVertexMultipleAttributeMatching(graph, parallelvertexidsets[for_core],  rulevertextuple.first);
            }

            for(int for_core=0;for_core<corenum;for_core++) {
                count+=matchedvertexids[for_core].size();
            }


            if(count>minsupport){

                Path temppath;
                temppath.pathidentifier=globalpathidcount;
                pathidentifier2pathsindex.insert({globalpathidcount,make_pair(0,paths[0].size())});
                globalpathidcount++;
                temppath.count=count;
                temppath.vertextuples.push_back(rulevertextuple.first);
                temppath.kleenestar=false;

                paths[0].push_back(temppath);

                #pragma omp parallel for
                for(int for_core=0;for_core<corenum;for_core++) {
                    for (auto matchedvertexid: matchedvertexids[for_core]) {

                        matchedvertices[matchedvertexid].pathidentifiers[0].push_back(temppath.pathidentifier);
                        vector<int> temptarget;
                        temptarget.push_back(matchedvertexid);
                        matchedvertices[matchedvertexid].targets[0].push_back(temptarget);
                    }
                }
                vector<int> allmatchedvertexids;
                allmatchedvertexids.clear();
                for(int for_core=0;for_core<corenum;for_core++) {
                    for(auto& vertexid: matchedvertexids[for_core]) {
                        allmatchedvertexids.push_back(vertexid);
                    }
                }
                matchedvertexs_attributeset.push_back(allmatchedvertexids);

                for(auto& pathidentifier: rulevertextuple.second){
                    paths[0][pathidentifier2pathsindex[pathidentifier].second].horizontalExtendingPathidentifiers.push_back(temppath.pathidentifier);
                }

//                paths[0][pathidentifier2pathsindex[path.pathidentifier].second].verticalExtendingPathidentifiers.push_back(temppath.pathidentifier);


            }
        }
        start=end;
        end=paths[0].size();
        numberOfCandidateAttribute++;
        if(start==end)break;
    }

    maximumsizemultipleattributes=paths[0].back().vertextuples[0].size();
}




void AssociationRules::ParallelFindFrequentLengthLPath(OriginalGraph& graph, int length,  vector<Path>& path_0, int maxattributesize) {


    int countmax=0;

    int test=0;

//    vector<vector<int>> parallelvertexidsets;
//    parallelvertexidsets.resize(corenum);
//    for(int i=0;i<graph.number_of_vertices;i++){
//        parallelvertexidsets[i%corenum].push_back(i);
//    }
    cout<<"START finding length "<< length << " frequent paths"<<endl;

    for(auto& path: path_0){//extends paths, we call it path A
        CheckTime();

//        bool continueFlag = false;
//        for (int i = 0; i < length; i++){
//            if (path.vertextuples[i].size() > maxattributesize) {
//                continueFlag = true;
//                break;
//            }
//        }
//        if(continueFlag)continue;

        //cout<<"test L: "<<test;
       // path.Show();

        int count = 0;
        vector<vector<int>> matchedvertexids;
        matchedvertexids.resize(corenum);

//        matchedvertexids.clear();

        vector<vector<pair<int,int>>> matchedvertexsets;
        matchedvertexsets.resize(corenum);

        #pragma omp parallel for
        for(int for_core=0;for_core<corenum;for_core++){
            matchedvertexsets[for_core].clear();
        }

        vector<int> pathid(corenum);
        #pragma omp parallel for
        for(int for_core=0;for_core<corenum;for_core++) {
            matchedvertexsets[for_core]=ParallelPathContainment(parallelvertexidsets[for_core], length-1, path.pathidentifier);
        }

//        int size=0;
//        for(int for_core=0;for_core<corenum;for_core++) {
//            size+=matchedvertexsets[for_core].size();
//        }
//         cout<<"test L: matchedsize="<<size<<endl;

        //int numberOfmatchedvertex = matchedvertex_pathidindexs.size();
        vector<bool> insertedtemptaregets;
        insertedtemptaregets.resize(graph.number_of_vertices);


        for(auto for_edgelabel: candidateedgelabels[length-1]) {//for all edge labels

            //cout<<for_edgelabel<<" / "<<candidateedgelabels[length-1].size()<<endl;
            vector<vector<int>> temptargets;
            vector<vector<vector<int>>> targetslist;
            targetslist.resize(corenum);

            #pragma omp parallel for
            for(int for_core=0;for_core<corenum;for_core++) {
                targetslist[for_core]=ParallelFindTarget(graph, matchedvertexsets[for_core], for_edgelabel, length);
            }

//            int size=0;
//            for(int for_core=0;for_core<corenum;for_core++) {
//                for(int i=0;i<targetslist[for_core].size();i++) {
//                    size += targetslist[for_core][i].size();
//                }
//            }
//            cout<<"test L: targetsize="<<size<<endl;

            vector<vector<int>> nexttargets;
            nexttargets.resize(graph.number_of_vertices);
            vector<bool> matchedvertexresult;
            matchedvertexresult.resize(graph.number_of_vertices);

            vector<bool> insertedmatchedvertexs;
            insertedmatchedvertexs.resize(graph.number_of_vertices);


            for(auto path0: targettuples_edges[length-1][for_edgelabel]) {//for all frequent vertex tuples
                if( path0.second.size() > maxattributesize)continue;

                test++;

                count=0;
                for(int i=0;i<graph.number_of_vertices;i++){
                    nexttargets[i].clear();
                    insertedmatchedvertexs[i]=false;
                    matchedvertexresult[i]=false;
                }

                #pragma omp parallel for
                for(int for_core=0;for_core<corenum;for_core++) {
                    matchedvertexids[for_core].clear();
                    matchedvertexids[for_core]=ParallelVertexMatching(graph, matchedvertexsets[for_core], targetslist[for_core], nexttargets, path0.first);
                }

                for(int for_core=0;for_core<corenum;for_core++) {
                    count+=matchedvertexids[for_core].size();
                }

                //cout<<"test L: "<<test<<"; edge="<<for_edgelabel<<", att="<<path0.second[0]<<" count="<<count<<endl;

                count*=1/ApproximateSamplingRate;

                if (count > minsupport) {

                    Path temppath;
                    temppath.pathidentifier = globalpathidcount;
                    pathidentifier2pathsindex.insert({globalpathidcount,make_pair(length,paths[length].size())});
                    globalpathidcount++;
                    temppath.count = count;
                    temppath.vertextuples = path.vertextuples;
                    temppath.vertextuples.push_back(path0.second);
                    temppath.edgelabels = path.edgelabels;
                    temppath.edgelabels.push_back(for_edgelabel);
                    temppath.kleenestar = false;

                    paths[length].push_back(temppath);

//                    cout<<"test frequentpath: ";
//                    temppath.Show();

                    #pragma omp parallel for
                    for(int for_core=0;for_core<corenum;for_core++) {
                        for (auto matchedvertexid: matchedvertexids[for_core]) {

                            matchedvertices[matchedvertexid].pathidentifiers[length].push_back(temppath.pathidentifier);
                            matchedvertices[matchedvertexid].targets[length].push_back(nexttargets[matchedvertexid]);
                        }
                    }


                    paths[length-1][pathidentifier2pathsindex[path.pathidentifier].second].verticalExtendingPathidentifiers.push_back(temppath.pathidentifier);
                }
            }
        }
    }
}



void AssociationRules::ParallelEfficientFindFrequentKleenePath(OriginalGraph& graph) {

    struct timespec temp1_startTime, temp1_endTime;
    struct timespec temp2_startTime, temp2_endTime;
    struct timespec temp3_startTime, temp3_endTime;
    struct timespec temp4_startTime, temp4_endTime;

    vector<bool> insertedtemptarget(graph.number_of_vertices,false);
    vector<bool> insertednexttarget(graph.number_of_vertices,false);
    vector<vector<int>> nexttargets;
    nexttargets.resize(graph.number_of_vertices);
    vector<bool> alltargets;
    alltargets.resize(graph.number_of_vertices);

//    vector<vector<int>> parallelvertexidsets;
//    parallelvertexidsets.resize(corenum);
//    for(int i=0;i<graph.number_of_vertices;i++){
//        parallelvertexidsets[i%corenum].push_back(i);
//    }

    int start= paths[1].size();
    int end;

    int doneedges=0;
    double testtime=0;

    //#pragma omp parallel for
    for (int i=0; i < candidateedgelabels[0].size();i++) {
        ParallelKleenePathMatching_perEdgeLabel(graph,candidateedgelabels[0][i]);
    }//for_edge

    end=paths[1].size();

    int numberOfCandidateAttribute=1;
    vector<pair<Path,vector<int>>>tempcandidatepaths;
    vector<pair<Path,vector<int>>>candidatepaths;
    //vector<vector<int>> nexttargets;
    nexttargets.resize(graph.number_of_vertices);


    //cout<<start<<","<<end<<endl;
    while(1) {

        tempcandidatepaths.clear();
        candidatepaths.clear();
        //Candidate generation

        for (int for_path1 = start; for_path1 < end - 1; for_path1++) {

            for (int for_path2 = for_path1 + 1; for_path2 < end; for_path2++) {
                Path candidatepath1;
                Path candidatepath2;
                Path candidatepath3;
                vector<int> pathidentifiers1;pathidentifiers1.clear();
                vector<int> pathidentifiers2;pathidentifiers2.clear();
                vector<int> pathidentifiers3;pathidentifiers3.clear();
                bool containFlag1=false;
                bool insertFlag1=true;
                bool containFlag2=false;
                bool insertFlag2=true;
                bool containFlag3=false;
                bool insertFlag4=true;


                candidatepath1 = MergeTwoOneLengthPathSource(paths[1][for_path1],paths[1][for_path2], numberOfCandidateAttribute);
                if(!candidatepath1.vertextuples.empty()) {

                    for (int for_att = 0; for_att < candidatepath1.vertextuples[0].size(); for_att++) {
                        containFlag1=false;
                        Path tmppath = candidatepath1;
                        tmppath.vertextuples[0].erase(tmppath.vertextuples[0].begin() + for_att);

                        for (int for_path = start; for_path < end; for_path++) {

                            if (paths[1][for_path] == tmppath) {
                                containFlag1 = true;
                                pathidentifiers1.push_back(paths[1][for_path].pathidentifier);
                                break;
                            }
                        }
                                                                                                                                                                               if (!containFlag1) {
                            insertFlag1 = false;
                            break;
                        }
                    }
                    if(insertFlag1)tempcandidatepaths.push_back(make_pair(candidatepath1,pathidentifiers1));
                }
                else insertFlag1=false;

                candidatepath2 = MergeTwoOneLengthPathTarget(paths[1][for_path1],paths[1][for_path2], numberOfCandidateAttribute);

                if(!candidatepath2.vertextuples.empty()) {
                    //candidatepath2.Show();
                    for (int for_att = 0; for_att < candidatepath2.vertextuples[1].size(); for_att++) {

                        containFlag2 = false;

                        Path tmppath = candidatepath2;
                        tmppath.vertextuples[1].erase(tmppath.vertextuples[1].begin() + for_att);
                        //tmppath.Show();
                        for (int for_path = start; for_path < end; for_path++) {

                            if (paths[1][for_path] == tmppath) {
                                //paths[1][for_path].Show();
                                pathidentifiers2.push_back(paths[1][for_path].pathidentifier);
                                containFlag2 = true;
                                break;
                            }
                        }
                        if (!containFlag2) {
                            insertFlag2 = false;
                            break;
                        }
                    }
                    if(insertFlag2)tempcandidatepaths.push_back(make_pair(candidatepath2,pathidentifiers2));
                }
                else insertFlag2=false;

                candidatepath3 = MergeTwoPath(paths[1][for_path1],paths[1][for_path2], numberOfCandidateAttribute);
                if(!candidatepath3.vertextuples.empty()&&insertFlag1&&insertFlag2){
                    pathidentifiers3=vectorUnion(pathidentifiers1,pathidentifiers2);
                    tempcandidatepaths.push_back(make_pair(candidatepath3,pathidentifiers3));
                }
            }
        }

        for(auto& tempcandidate: tempcandidatepaths){
            //tempcandidate.first.Show();
            bool duplicateFlag=false;
            for(auto& candidate: candidatepaths){

                if(candidate.first == tempcandidate.first){
                    duplicateFlag=true;
                    break;
                }
            }
            if(!duplicateFlag){
                //tempcandidate.first.Show();
                candidatepaths.push_back(tempcandidate);
            }

        }

        //cout<<"kleene candi: "<<candidatepaths.size()<<endl;

        //Find frequent Set//
        int test=0;
        for(auto& rulepath: candidatepaths){
            test++;
            //cout<<"test"<<test<<endl;

            int count = 0;

            vector<vector<int>> matchedvertexsets;
            matchedvertexsets.resize(corenum);

            count=0;
            for(int i=0;i<graph.number_of_vertices;i++){
                nexttargets[i].clear();
            }

            #pragma omp parallel for
            for(int for_core=0;for_core<corenum;for_core++){
                matchedvertexsets[for_core]=ParallelVertexMatching_forEfficient(graph, parallelvertexidsets[for_core], rulepath.second,1, nexttargets);
            }

            for(int for_core=0;for_core<corenum;for_core++){
                count+=matchedvertexsets[for_core].size();
            }
            count*=1/ApproximateSamplingRate;

            if(count>minsupport){

                Path temppath;
                temppath.pathidentifier=globalpathidcount;
                pathidentifier2pathsindex.insert({globalpathidcount,make_pair(1,paths[1].size())});
                globalpathidcount++;
                temppath.count=count;
                temppath.vertextuples=rulepath.first.vertextuples;
                temppath.edgelabels=rulepath.first.edgelabels;
                temppath.kleenestar=true;

                paths[1].push_back(temppath);

                #pragma omp parallel for
                for(int for_core=0;for_core<corenum;for_core++) {
                    for (auto matchedvertexid: matchedvertexsets[for_core]) {

                        matchedvertices[matchedvertexid].pathidentifiers[1].push_back(temppath.pathidentifier);
                        matchedvertices[matchedvertexid].targets[1].push_back(nexttargets[matchedvertexid]);
                    }
                }
                for(auto& pathidentifier: rulepath.second){
                    paths[1][pathidentifier2pathsindex[pathidentifier].second].horizontalExtendingPathidentifiers.push_back(temppath.pathidentifier);
                }
            }
        }
        start=end;
        end=paths[1].size();
        numberOfCandidateAttribute++;
        if(start==end)break;
    }
}


void AssociationRules::ParallelKleenePathMatching_perEdgeLabel(OriginalGraph& graph, int edgelabel){


    struct timespec temp1_startTime, temp1_endTime;
    struct timespec temp2_startTime, temp2_endTime;
    struct timespec temp3_startTime, temp3_endTime;
    struct timespec temp4_startTime, temp4_endTime;

    vector<bool> insertedtemptarget(graph.number_of_vertices,false);
    vector<bool> insertednexttarget(graph.number_of_vertices,false);
    vector<vector<int>> nexttargets;
    nexttargets.resize(graph.number_of_vertices);
    vector<bool> alltargets;
    alltargets.resize(graph.number_of_vertices);
    clock_gettime(CLOCK_REALTIME, &temp2_startTime);
    vector<vector<int>> edgetargetlists;
    edgetargetlists.clear();
    vector<bool> bfsdoneflag;

    for (int i = 0; i < graph.number_of_vertices; i++) {
        alltargets[i] = false;
    }

    edgetargetlists.resize(graph.number_of_vertices);
    bfsdoneflag.resize(graph.number_of_vertices);
    #pragma omp parallel for
    for (int for_core = 0; for_core < corenum; for_core++) {
        for (auto &vertexid : parallelvertexidsets[for_core]) {
            edgetargetlists[vertexid].clear();
            bfsdoneflag[vertexid]=false;
        }
    }

    vector<int> temptargets;
    vector<int> bfssearchvertex;
    vector<int> bfssearchvertex_next;
    int parallel_for_length;
    int parallel_i;
    int parallel_edgeid;
    Edge tempedge;
    int for_core_vertex;
    int bfsvertexid;
    int for_target;
    int parallel_vertexid;

    //cout<<"test size: "<<insertedtemptarget.size()<<endl;

    #pragma omp parallel for private(insertedtemptarget, temptargets, bfssearchvertex, bfssearchvertex_next, parallel_for_length, parallel_i, parallel_edgeid, tempedge, for_core_vertex, bfsvertexid, for_target,parallel_vertexid)
    for (int for_core = 0; for_core < corenum; for_core++) {
        for (for_core_vertex = 0; for_core_vertex < parallelvertexidsets[for_core].size(); for_core_vertex++) {

            bfsvertexid = parallelvertexidsets[for_core][for_core_vertex];
            //cout<<"bfs for_vertex1 : "<<bfsvertexid<<endl;
            temptargets.clear();
            bfssearchvertex.clear();

            // the vertex is unnecessary for check//
            if (matchedvertices[bfsvertexid].pathidentifiers[0].empty() || !graph.vertex_connectingcount_edgelabels[bfsvertexid][edgelabel]) {
                continue;
            }

            bfssearchvertex.push_back(bfsvertexid);

            insertedtemptarget.resize(graph.number_of_vertices);
            for (parallel_i = 0; parallel_i < graph.number_of_vertices; parallel_i++) {
                //insertednexttarget[i]=false;
                //cout<<"paralllel i = "<<insertedtemptarget.size()<<","<<parallel_i<<endl;
                insertedtemptarget[parallel_i] = false;
            }
            insertedtemptarget[bfsvertexid] = true;
            bfssearchvertex_next.clear();
            //cout<<"bfs for_vertex2 : "<<bfsvertexid<<endl;
            for (parallel_for_length = 0; parallel_for_length < maximumpathlength; parallel_for_length++) {
                //bfssearchvertex_next.clear();
                while (1) {
                    if (bfssearchvertex.empty())break;
                    parallel_vertexid = bfssearchvertex.back();
                    bfssearchvertex.pop_back();


                    for (parallel_edgeid=0; parallel_edgeid < graph.vertex_connecting_edge_list[parallel_vertexid].size();parallel_edgeid++) {

                        tempedge = graph.edge_list[graph.vertex_connecting_edge_list[parallel_vertexid][parallel_edgeid]];
                        if (tempedge.edgelabel == edgelabel) {

                            if (!insertedtemptarget[tempedge.dst]) {

                                bfssearchvertex_next.push_back(tempedge.dst);
                                temptargets.push_back(tempedge.dst);
                                insertedtemptarget[tempedge.dst] = true;
                                alltargets[tempedge.dst] = true;
                            }
                        }
                    }
                }
                if (bfssearchvertex_next.empty())break;
                bfssearchvertex = bfssearchvertex_next;
                bfssearchvertex_next.clear();
            }
            edgetargetlists[bfsvertexid]=temptargets;
            bfsdoneflag[bfsvertexid]=true;

        }
    }

    clock_gettime(CLOCK_REALTIME, &temp2_endTime);
    cout<<"---BFS :"<<edgelabel<<"), time= ";
    cout<<(temp2_endTime.tv_sec - temp2_startTime.tv_sec) + (temp2_endTime.tv_nsec - temp2_startTime.tv_nsec) * pow(10, -9)<<endl;

    vector<bool> insertedmatchedvertexs;
    vector<bool> matchedvertexresult;
    matchedvertexresult.resize(graph.number_of_vertices);
    insertedmatchedvertexs.resize(graph.number_of_vertices);

    int testcount=0;
    for (auto& sourceattribute: paths[0]) {// for the start of paths
        if(sourceattribute.vertextuples[0].size()>1)continue;

        vector<int> matchedsourceids;
        matchedsourceids.clear();

        vector<int> pathmatchedvertexids;

        for (auto targetattributes: targettuples_edges[1][edgelabel]) {// for the end of paths
            if(targetattributes.second.size() >1)continue;

            int count = 0;
            for (int i = 0; i < graph.number_of_vertices; i++) {
                nexttargets[i].clear();
                insertedmatchedvertexs[i] = false;
                matchedvertexresult[i]=false;
            }
            pathmatchedvertexids.clear();
            std::mutex mtx;
            int matchedvertexid;
            int for_target;

            //cout<<"testcount: "<<testcount<<endl;
            #pragma omp parallel for private(matchedvertexid, for_target)
            for (int for_vertex=0; for_vertex<matchedvertexs_attributeset[sourceattribute.pathidentifier].size();for_vertex++) {
                matchedvertexid =  matchedvertexs_attributeset[sourceattribute.pathidentifier][for_vertex];

                for (for_target=0; for_target<edgetargetlists[matchedvertexid].size();for_target++) {

                    if (Contains(matchedvertices[edgetargetlists[matchedvertexid][for_target]].targetidentifiers, targetattributes.first)) {

                        if (!insertedmatchedvertexs[matchedvertexid]) {
                            mtx.lock();
                            pathmatchedvertexids.push_back(matchedvertexid);
                            count++;
                            mtx.unlock();
                            insertedmatchedvertexs[matchedvertexid] = true;
                        }
                        nexttargets[matchedvertexid].push_back(edgetargetlists[matchedvertexid][for_target]);
                    }
                }
            }
            //cout<<"count : "<<count<<endl;

            count*=1/ApproximateSamplingRate;

            if (count > minsupport) {

                Path temppath;
                temppath.pathidentifier = globalpathidcount;
                temppath.count = count;
                temppath.vertextuples.push_back(sourceattribute.vertextuples[0]);
                temppath.vertextuples.push_back(targetattributes.second);
                temppath.edgelabels.push_back(edgelabel);
                temppath.kleenestar = true;

                mtx.lock();
                pathidentifier2pathsindex.insert({globalpathidcount,make_pair(1,paths[1].size())});
                globalpathidcount++;
                paths[1].push_back(temppath);

                for (auto matchedvertexid: pathmatchedvertexids) {
                    matchedvertices[matchedvertexid].pathidentifiers[1].push_back(temppath.pathidentifier);
                    matchedvertices[matchedvertexid].targets[1].push_back(nexttargets[matchedvertexid]);
                }

                sourceattribute.verticalExtendingPathidentifiers.push_back(temppath.pathidentifier);
                mtx.unlock();

                //testcount++;
            }

        }//for_target candidate

    }//for_source

}



void AssociationRules::ParallelFindFrequentKleenePath(OriginalGraph& graph) {


    struct timespec temp1_startTime, temp1_endTime;
    struct timespec temp2_startTime, temp2_endTime;
    struct timespec temp3_startTime, temp3_endTime;
    struct timespec temp4_startTime, temp4_endTime;

    vector<bool> insertedtemptarget(graph.number_of_vertices,false);
    vector<bool> insertednexttarget(graph.number_of_vertices,false);
    vector<vector<int>> nexttargets;
    nexttargets.resize(graph.number_of_vertices);
    vector<bool> alltargets;
    alltargets.resize(graph.number_of_vertices);
    clock_gettime(CLOCK_REALTIME, &temp2_startTime);
    vector<vector<int>> edgetargetlists;
    edgetargetlists.clear();
    vector<bool> bfsdoneflag;

    for (int i = 0; i < graph.number_of_vertices; i++) {
        alltargets[i] = false;
    }


    edgetargetlists.resize(graph.number_of_vertices);
    bfsdoneflag.resize(graph.number_of_vertices);
//    #pragma omp parallel for
//    for (int for_core = 0; for_core < corenum; for_core++) {
//        for (auto &for_vertex : parallelvertexidsets[for_core]) {
//            edgetargetlists[for_vertex].clear();
//            bfsdoneflag[for_vertex]=false;
//        }
//    }

    vector<int> temptargets;
    vector<int> bfssearchvertex;
    vector<int> bfssearchvertex_next;
    int parallel_for_length;
    int parallel_i;
    int parallel_edgeid;
    Edge tempedge;
    int for_core_vertex;
    int bfsvertexid;
    int for_target;
    int parallel_vertexid;


    int doneedges=0;
    double testtime=0;
    for (auto for_edgelabel: candidateedgelabels[0]) {

        clock_gettime(CLOCK_REALTIME, &temp2_startTime);
//        vector<vector<int>> edgetargetlists;
//        edgetargetlists.clear();
        #pragma omp parallel for
        for (int for_core = 0; for_core < corenum; for_core++) {
            for (auto &for_vertex : parallelvertexidsets[for_core]) {
                edgetargetlists[for_vertex].clear();
                bfsdoneflag[for_vertex]=false;
            }
        }

        for (int i = 0; i < graph.number_of_vertices; i++) {
            alltargets[i] = false;
        }
        #pragma omp parallel for private(insertedtemptarget, temptargets, bfssearchvertex, bfssearchvertex_next, parallel_for_length, parallel_i, parallel_edgeid, tempedge, for_core_vertex, bfsvertexid, for_target,parallel_vertexid)
        for (int for_core = 0; for_core < corenum; for_core++) {
            for (for_core_vertex = 0; for_core_vertex < parallelvertexidsets[for_core].size(); for_core_vertex++) {
                //cout<<corenum<<", test4"<<endl;
                //for (int bfsvertexid = 0; bfsvertexid < graph.number_of_vertices; bfsvertexid++) {

                bfsvertexid = parallelvertexidsets[for_core][for_core_vertex];
                //cout << "start: " << for_core << "," << bfsvertexid << endl;
                vector<int> temptargets;
                vector<int> bfssearchvertex;
                temptargets.clear();
                bfssearchvertex.clear();

                // the vertex is unnecessary for check//
                if (matchedvertices[bfsvertexid].pathidentifiers[0].empty() || !graph.vertex_connectingcount_edgelabels[bfsvertexid][for_edgelabel]) {
                    //edgetargetlists.push_back(temptargets);
                    continue;
                }

                bfssearchvertex.push_back(bfsvertexid);

                insertedtemptarget.resize(graph.number_of_vertices);
                for (parallel_i = 0; parallel_i < graph.number_of_vertices; parallel_i++) {
                    //insertednexttarget[i]=false;
                    //cout<<"paralllel i = "<<insertedtemptarget.size()<<","<<parallel_i<<endl;
                    insertedtemptarget[parallel_i] = false;
                }

                insertedtemptarget[bfsvertexid] = true;
                bfssearchvertex_next.clear();

                for (parallel_for_length = 0; parallel_for_length < maximumpathlength; parallel_for_length++) {
                    //bfssearchvertex_next.clear();
                    while (1) {
                        if (bfssearchvertex.empty())break;
                        parallel_vertexid = bfssearchvertex.back();
                        bfssearchvertex.pop_back();

                        for (parallel_edgeid=0; parallel_edgeid < graph.vertex_connecting_edge_list[parallel_vertexid].size();parallel_edgeid++) {

                            tempedge = graph.edge_list[graph.vertex_connecting_edge_list[parallel_vertexid][parallel_edgeid]];
                            if (tempedge.edgelabel == for_edgelabel) {

                                if (!insertedtemptarget[tempedge.dst]) {
                                    //if (tempedge.dst < bfsvertexid) {
//                                        temptargets.push_back(tempedge.dst);
//                                        insertedtemptarget[tempedge.dst] = true;
//                                        for (for_target=0; for_target < edgetargetlists[tempedge.dst].size();for_target++) {
//                                            if (!insertedtemptarget[edgetargetlists[tempedge.dst][for_target]]) {
//                                            temptargets.push_back(edgetargetlists[tempedge.dst][for_target]);
//                                            insertedtemptarget[edgetargetlists[tempedge.dst][for_target]] = true;
//                                        }
//                                     }
//                                    } else {
                                        bfssearchvertex_next.push_back(tempedge.dst);
                                        temptargets.push_back(tempedge.dst);
                                        insertedtemptarget[tempedge.dst] = true;
                                        alltargets[tempedge.dst] = true;
                                    //}
                                }
                            }
                        }
                    }
                    if (bfssearchvertex_next.empty())break;
                    bfssearchvertex = bfssearchvertex_next;
                    bfssearchvertex_next.clear();
                }
                edgetargetlists[bfsvertexid]=temptargets;
                bfsdoneflag[bfsvertexid]=true;

                //clock_gettime(CLOCK_REALTIME, &frequentkleenepath_endTime);
                //cout<<"---DURING BFS ("<<bfsvertexid<<"/"<< graph.number_of_vertices<<", size="<<temptargets.size()<<"): ";
                //cout<<(frequentkleenepath_endTime.tv_sec - frequentkleenepath_startTime.tv_sec) + (frequentkleenepath_endTime.tv_nsec - frequentkleenepath_startTime.tv_nsec) * pow(10, -9)<<endl;
            }
        }

        doneedges++;
        clock_gettime(CLOCK_REALTIME, &temp2_endTime);
        cout<<"---BFS :"<<for_edgelabel<<" ("<<doneedges<<"/"<< candidateedgelabels[0].size()<<"), time= ";
        cout<<(temp2_endTime.tv_sec - temp2_startTime.tv_sec) + (temp2_endTime.tv_nsec - temp2_startTime.tv_nsec) * pow(10, -9)<<endl;

        //clock_gettime(CLOCK_REALTIME, &temp1_startTime);

        vector<bool> insertedmatchedvertexs;
        vector<bool> matchedvertexresult;
        matchedvertexresult.resize(graph.number_of_vertices);
        insertedmatchedvertexs.resize(graph.number_of_vertices);

        int testcount=0;
        for (auto& sourceattributes: paths[0]) {// for the start of paths

            vector<int> matchedsourceids;
            matchedsourceids.clear();

            vector<int> pathmatchedvertexids;

            for (auto& targetattributes: targettuples_edges[1][for_edgelabel]) {// for the end of paths
                int count = 0;
                for (int i = 0; i < graph.number_of_vertices; i++) {
                    nexttargets[i].clear();
                    insertedmatchedvertexs[i] = false;
                    matchedvertexresult[i]=false;
                }
                pathmatchedvertexids.clear();
                std::mutex mtx;
                int matchedvertexid;
                int for_target;
                vector<int> matchedsourceids;
                matchedsourceids.clear();


                #pragma omp parallel for private(matchedvertexid, for_target)
                for (int for_vertex=0; for_vertex<matchedvertexs_attributeset[sourceattributes.pathidentifier].size();for_vertex++) {
                    matchedvertexid =  matchedvertexs_attributeset[sourceattributes.pathidentifier][for_vertex];

                    for (for_target=0; for_target<edgetargetlists[matchedvertexid].size();for_target++) {

                        if (Contains(matchedvertices[edgetargetlists[matchedvertexid][for_target]].targetidentifiers, targetattributes.first)) {

                            if (!insertedmatchedvertexs[matchedvertexid]) {
                                mtx.lock();
                                pathmatchedvertexids.push_back(matchedvertexid);
                                count++;
                                mtx.unlock();
                                insertedmatchedvertexs[matchedvertexid] = true;
                            }
                            nexttargets[matchedvertexid].push_back(edgetargetlists[matchedvertexid][for_target]);
                        }
                    }
                }

                count*=1/ApproximateSamplingRate;

                if (count > minsupport) {

                    Path temppath;
                    temppath.pathidentifier = globalpathidcount;

                    temppath.count = count;
                    temppath.vertextuples.push_back(sourceattributes.vertextuples[0]);
                    temppath.vertextuples.push_back(targetattributes.second);
                    temppath.edgelabels.push_back(for_edgelabel);
                    temppath.kleenestar = true;

                    mtx.lock();
                    pathidentifier2pathsindex.insert({globalpathidcount,make_pair(1,paths[1].size())});
                    globalpathidcount++;
                    paths[1].push_back(temppath);

                    for (auto matchedvertexid: pathmatchedvertexids) {
                        matchedvertices[matchedvertexid].pathidentifiers[1].push_back(temppath.pathidentifier);
                        matchedvertices[matchedvertexid].targets[1].push_back(nexttargets[matchedvertexid]);
                    }
                    mtx.unlock();
                    testcount++;
                }

            }//for_target candidate

        }//for_source
        //cout<<testtime<<endl;
        //clock_gettime(CLOCK_REALTIME, &temp1_endTime);
        //cout<<"---Matching ("<<doneedges<<"/"<< candidateedgelabels.size()<<"), time= ";
        //cout<<(temp1_endTime.tv_sec - temp1_startTime.tv_sec) + (temp1_endTime.tv_nsec - temp1_startTime.tv_nsec) * pow(10, -9)<<endl;

        //doneedges++;
        //clock_gettime(CLOCK_REALTIME, &frequentkleenepath_endTime);
        //cout<<"During ("<<doneedges<<"/"<< candidateedgelabels.size()<<"), time= ";
        //cout<<(frequentkleenepath_endTime.tv_sec - frequentkleenepath_startTime.tv_sec) + (frequentkleenepath_endTime.tv_nsec - frequentkleenepath_startTime.tv_nsec) * pow(10, -9)<<endl;

    }//for_edge

}


//void AssociationRules::FindCandidateEdgeLabel(OriginalGraph& graph) {
//
//    candidateedgelabels.clear();
//
//    for(int i=0;i<graph.cardinality_of_edgelabels;i++){
//        if(graph.numbers_of_edgelabels[i]>minsupport)candidateedgelabels.push_back(i);
//    }
//    cout<<"candidate edge size = "<<candidateedgelabels.size()<<endl;
//}




void AssociationRules::ParallelEfficientFindFrequentLengthOnePath(OriginalGraph &graph){

    struct timespec temp1_startTime, temp1_endTime;
    double temp1time;
    clock_gettime(CLOCK_REALTIME, &temp1_startTime);

    vector<Path> paths_0;
    //paths_0.resize(maximumsizemultipleattributes);

    cout<<"max attributesize: "<<maximumsizemultipleattributes<<endl;

    //for(int i=0;i<maximumsizemultipleattributes;i++){
    paths_0.clear();
    //}
    for(int i=0;i<paths[0].size();i++){
        if(EmptysourceFlag&&paths[0][i].vertextuples[0].size()==0)paths_0.push_back(paths[0][i]);
        if(paths[0][i].vertextuples[0].size()==1)paths_0.push_back(paths[0][i]);
    }


    ParallelFindFrequentLengthOnePath(graph,paths_0,1);

    clock_gettime(CLOCK_REALTIME, &temp1_endTime);
    temp1time=(temp1_endTime.tv_sec - temp1_startTime.tv_sec) + (temp1_endTime.tv_nsec - temp1_startTime.tv_nsec) * pow(10, -9);
    cout << "unit length 1 paths: "<< "size" <<paths[1].size()<<", time "<<temp1time<<endl;

    int start=0;
    int end=paths[1].size();
    int pathlength=1;
    int maximumsizecombination=paths[1].size();
    int numberOfCandidateAttribute=1;
    vector<pair<Path,vector<int>>>tempcandidatepaths;
    vector<pair<Path,vector<int>>>candidatepaths;
    vector<vector<int>> nexttargets;
    nexttargets.resize(graph.number_of_vertices);

//    vector<vector<int>> parallelvertexidsets;
//    parallelvertexidsets.resize(corenum);
//    for(int i=0;i<graph.number_of_vertices;i++){
//        parallelvertexidsets[i%corenum].push_back(i);
//    }

    while(1) {

        tempcandidatepaths.clear();
        candidatepaths.clear();
        //Candidate generation

        for (int for_path1 = start; for_path1 < end - 1; for_path1++) {

            for (int for_path2 = for_path1 + 1; for_path2 < end; for_path2++) {
                Path candidatepath1; candidatepath1.clear();
                Path candidatepath2; candidatepath2.clear();
                Path candidatepath3; candidatepath3.clear();
                vector<int> pathidentifiers1; pathidentifiers1.clear();
                vector<int> pathidentifiers2; pathidentifiers2.clear();
                vector<int> pathidentifiers3; pathidentifiers3.clear();
                bool containFlag1=false;
                bool insertFlag1=true;
                bool containFlag2=false;
                bool insertFlag2=true;
                bool containFlag3=false;
                bool insertFlag4=true;

                if(!EmptysourceFlag&&(paths[1][for_path1].vertextuples[0].empty()||paths[1][for_path2].vertextuples[0].empty())){
                    candidatepath1 = MergeTwoOneLengthPathSource(paths[1][for_path1],paths[1][for_path2], numberOfCandidateAttribute);
                }
                if(!candidatepath1.vertextuples.empty()) {

                    for (int for_att = 0; for_att < candidatepath1.vertextuples[0].size(); for_att++) {
                        containFlag1=false;
                        Path tmppath = candidatepath1;
                        tmppath.vertextuples[0].erase(tmppath.vertextuples[0].begin() + for_att);

                        for (int for_path = start; for_path < end; for_path++) {

                            if (paths[1][for_path] == tmppath) {
                                containFlag1 = true;
                                pathidentifiers1.push_back(paths[1][for_path].pathidentifier);
                                break;
                            }
                        }
                        if (!containFlag1) {
                            insertFlag1 = false;
                            break;
                        }
                    }
                    if(insertFlag1)tempcandidatepaths.push_back(make_pair(candidatepath1,pathidentifiers1));
                }
                else insertFlag1=false;

                candidatepath2 = MergeTwoOneLengthPathTarget(paths[1][for_path1],paths[1][for_path2], numberOfCandidateAttribute);

                if(!candidatepath2.vertextuples.empty()) {
                    //candidatepath2.Show();
                    for (int for_att = 0; for_att < candidatepath2.vertextuples[1].size(); for_att++) {

                        containFlag2 = false;

                        Path tmppath = candidatepath2;
                        tmppath.vertextuples[1].erase(tmppath.vertextuples[1].begin() + for_att);
                        //tmppath.Show();
                        for (int for_path = start; for_path < end; for_path++) {

                            if (paths[1][for_path] == tmppath) {
                                //paths[1][for_path].Show();
                                pathidentifiers2.push_back(paths[1][for_path].pathidentifier);
                                containFlag2 = true;
                                break;
                            }
                        }
                        if (!containFlag2) {
                            insertFlag2 = false;
                            break;
                        }
                    }
                    if(insertFlag2)tempcandidatepaths.push_back(make_pair(candidatepath2,pathidentifiers2));
                }
                else insertFlag2=false;

                if(!EmptysourceFlag&&(paths[1][for_path1].vertextuples[0].empty()||paths[1][for_path2].vertextuples[0].empty())) {
                    candidatepath3 = MergeTwoPath(paths[1][for_path1], paths[1][for_path2], numberOfCandidateAttribute);
                }
                if(!candidatepath3.vertextuples.empty()&&insertFlag1&&insertFlag2){
                    pathidentifiers3=vectorUnion(pathidentifiers1,pathidentifiers2);
                    tempcandidatepaths.push_back(make_pair(candidatepath3,pathidentifiers3));
                }
            }
        }

        for(auto& tempcandidate: tempcandidatepaths){
            //tempcandidate.first.Show();
            bool duplicateFlag=false;
            for(auto& candidate: candidatepaths){

                if(candidate.first == tempcandidate.first){
                    duplicateFlag=true;
                    break;
                }
            }
            if(!duplicateFlag){
                //tempcandidate.first.Show();
                candidatepaths.push_back(tempcandidate);
            }

        }

        clock_gettime(CLOCK_REALTIME, &temp1_endTime);
        temp1time=(temp1_endTime.tv_sec - temp1_startTime.tv_sec) + (temp1_endTime.tv_nsec - temp1_startTime.tv_nsec) * pow(10, -9);
        //cout << "temp1 = "<<temp1time<<endl;

        //Find frequent Set//
        int test=0;
        cout<<"candidate size: "<< candidatepaths.size()<<endl;

        for(auto& rulepath: candidatepaths){
            test++;
            //cout<<"test"<<test<<endl;

            int count = 0;


            count=0;
            for(int i=0;i<graph.number_of_vertices;i++){
                nexttargets[i].clear();
            }



            vector<vector<int>> matchedvertexsets;
            matchedvertexsets.resize(corenum);

            #pragma omp parallel for
            for(int for_core=0;for_core<corenum;for_core++){
                matchedvertexsets[for_core]=ParallelVertexMatching_forEfficient(graph, parallelvertexidsets[for_core], rulepath.second,1, nexttargets);
            }

            for(int for_core=0;for_core<corenum;for_core++){
                count+=matchedvertexsets[for_core].size();
            }

            //cout<<"test count:= " <<count<<endl;
            count*=1/ApproximateSamplingRate;
            //cout<<"test aproximate count:= "<< count<<endl;

            if(count>minsupport){

                Path temppath;
                temppath.pathidentifier=globalpathidcount;
                pathidentifier2pathsindex.insert({globalpathidcount,make_pair(1,paths[1].size())});
                globalpathidcount++;
                temppath.count=count;
                temppath.vertextuples=rulepath.first.vertextuples;
                temppath.edgelabels=rulepath.first.edgelabels;
                temppath.kleenestar=false;

                paths[1].push_back(temppath);

                #pragma omp parallel for
                for(int for_core=0;for_core<corenum;for_core++) {
                    for (auto matchedvertexid: matchedvertexsets[for_core]) {

                        matchedvertices[matchedvertexid].pathidentifiers[1].push_back(temppath.pathidentifier);
                        matchedvertices[matchedvertexid].targets[1].push_back(nexttargets[matchedvertexid]);
                    }
                }

                for(auto& pathidentifier: rulepath.second){
                    paths[1][pathidentifier2pathsindex[pathidentifier].second].horizontalExtendingPathidentifiers.push_back(temppath.pathidentifier);
                }

            }
        }
        start=end;
        end=paths[1].size();
        numberOfCandidateAttribute++;
        if(start==end)break;
    }

}

void AssociationRules::ParallelEfficientFindFrequentLengthLPath(OriginalGraph &graph, int length){

    struct timespec temp1_startTime, temp1_endTime;
    double temp1time;
    clock_gettime(CLOCK_REALTIME, &temp1_startTime);

    vector<Path> paths_0;
    paths_0.clear();
    //paths_0.resize(maximumsizemultipleattributes);

    cout<<"max attributesize: "<<maximumsizemultipleattributes<<endl;

//    for(int i=0;i<maximumsizemultipleattributes;i++){
//        paths_0.clear();
//    }
    for(int i=0;i<paths[length-1].size();i++){
        bool insertFlag=true;
        for(int j=0;j<length;j++){
            if(paths[length-1][i].vertextuples[j].size()>1){
                insertFlag=false;
                break;
            }
//            if(j==length-1&&graph.attribute_connecting_vertex[paths[length-1][i].vertextuples[length-1][0]].size()<=minsupport){
//                insertFlag=false;
//            }
        }
        if(insertFlag)paths_0.push_back(paths[length-1][i]);
        //if(paths[length-1][i].vertextuples[0].size()==1)paths_0.push_back(paths[length-1][i]);
    }

    //cout<<"test2"<<endl;
    ParallelFindFrequentLengthLPath(graph,length,paths_0,1);
    //cout<<"test2"<<endl;

    clock_gettime(CLOCK_REALTIME, &temp1_endTime);
    temp1time=(temp1_endTime.tv_sec - temp1_startTime.tv_sec) + (temp1_endTime.tv_nsec - temp1_startTime.tv_nsec) * pow(10, -9);
    cout << "temp1 = "<<temp1time<<endl;

    int pathlength=1;
    int maximumsizecombination=paths[1].size();
    int numberOfCandidateAttribute=1;
    vector<pair<Path,vector<int>>>tempcandidatepaths;
    vector<pair<Path,vector<int>>>candidatepaths;
    vector<vector<int>> nexttargets;
    nexttargets.resize(graph.number_of_vertices);

//    vector<vector<int>> parallelvertexidsets;
//    parallelvertexidsets.resize(corenum);
//    for(int i=0;i<graph.number_of_vertices;i++){
//        parallelvertexidsets[i%corenum].push_back(i);
//    }

    int start=0;
    int end=paths[length].size();

    while(1) {

        tempcandidatepaths.clear();
        candidatepaths.clear();
        std::mutex mtx;
        //Candidate generation

        int for_path2;
        Path candidatepath;
        vector<int> pathidentifiers;
        vector<vector<int> > combinations;
        vector<int> enumerate;
        int i;
        int for_edgelabel;
        int for_length;
        Path path1;
        Path path2;
        vector<int> mergetuples;
        bool mergePosibilityFlag;

        #pragma omp parallel for private(for_path2,candidatepath,pathidentifiers,combinations,enumerate,i,for_edgelabel,for_length,path1,path2,mergePosibilityFlag,mergetuples)
        for (int for_path1 = start; for_path1 < end - 1; for_path1++) {

            //cout << for_path1 << "/" << end << endl;
            for (for_path2 = for_path1 + 1; for_path2 < end; for_path2++) {


                pathidentifiers.clear();

                //computecombination

                enumerate.clear();
                for (i = 0; i < length + 1; i++) {
                    enumerate.push_back(i);
                }
                //ShowVector(enumerate);

                //paths[length][for_path1].Show();
                //paths[length][for_path2].Show();

                //Compute Merge posibility;
                path1 = paths[length][for_path1];
                path2 = paths[length][for_path2];

                mergePosibilityFlag = true;
                if (path1.kleenestar != path2.kleenestar)mergePosibilityFlag = false;
                if (path1.edgelabels.size() != path2.edgelabels.size())mergePosibilityFlag = false;
                else {
                    for (for_edgelabel = 0; for_edgelabel < path1.edgelabels.size(); for_edgelabel++) {
                        if (path1.edgelabels[for_edgelabel] != path2.edgelabels[for_edgelabel]) {
                            mergePosibilityFlag = false;
                            break;
                        }
                    }
                }

                if (!mergePosibilityFlag)continue;

                for (for_length = 0; for_length < length + 1; for_length++) {
                    //compute combination (for_length)
                    combination(enumerate, for_length + 1, combinations);

//                    for(int i=0;i<combinations.size();i++){
//                        for(int j=0; j< combinations[i].size();j++){
//                            cout<<combinations[i][j]<< " ";
//                        }
//                        cout<<endl;
//                    }

                    while (!combinations.empty()) {


                        mergetuples = combinations.back();
                        combinations.pop_back();

                        //ShowVector(mergetuples);

                        candidatepath = MergeTwoPathGivenTuples(paths[length][for_path1], paths[length][for_path2], numberOfCandidateAttribute, mergetuples);

                        //candidatepath.Show();

                        if (!candidatepath.vertextuples.empty()) {
                            pathidentifiers.clear();
                            CandidatePathMatching(candidatepath, pathidentifiers, mergetuples, 0, start, end, length);

                            //ShowVector(pathidentifiers);
                            if (!pathidentifiers.empty()){
                                mtx.lock();
                                tempcandidatepaths.push_back(make_pair(candidatepath, pathidentifiers));
                                mtx.unlock();
                            }


                        }
                    }
                }
            }
        }


        cout<<"temp candidate size = "<<tempcandidatepaths.size()<<endl;

//        #pragma omp parallel shared(candidatepaths)
//        {
        bool duplicateFlag = false;
        int j;
        #pragma omp parallel for private(duplicateFlag,j)
        for (int i = 0; i < tempcandidatepaths.size(); i++) {
            //tempcandidate.first.Show();
            duplicateFlag = false;
            //cout << i <<"/"<<tempcandidatepaths.size()<<endl;
            for (j = i+1; j < tempcandidatepaths.size();j++) {

                if (tempcandidatepaths[j].first == tempcandidatepaths[i].first) {
                    duplicateFlag = true;
                    break;
                }
            }
            if (!duplicateFlag) {
                //tempcandidate.first.Show();
                mtx.lock();
                candidatepaths.push_back(tempcandidatepaths[i]);
                mtx.unlock();
            }
        }
        //}



//        for(auto& tempcandidate: tempcandidatepaths){
//            //tempcandidate.first.Show();
//            bool duplicateFlag=false;
//            for(auto& candidate: candidatepaths){
//
//                if(candidate.first == tempcandidate.first){
//                    duplicateFlag=true;
//                    break;
//                }
//            }
//            if(!duplicateFlag){
//                //tempcandidate.first.Show();
//                candidatepaths.push_back(tempcandidate);
//            }
//        }

        cout<<length<<" path candidate size: "<<candidatepaths.size()<<endl;

        //return;
        //if(candidatepaths.size())break;

        //Find frequent Set//
        CheckTime();
        int test=0;
        for(auto& rulepath: candidatepaths){
            test++;


            int count = 0;

            vector<vector<int>> matchedvertexsets;
            matchedvertexsets.resize(corenum);

            count=0;
            for(int i=0;i<graph.number_of_vertices;i++){
                nexttargets[i].clear();
            }

            #pragma omp parallel for
            for(int for_core=0;for_core<corenum;for_core++){
                matchedvertexsets[for_core]=ParallelVertexMatching_forEfficient(graph, parallelvertexidsets[for_core], rulepath.second,length, nexttargets);
            }

            for(int for_core=0;for_core<corenum;for_core++){
                count+=matchedvertexsets[for_core].size();
            }
            count*=1/ApproximateSamplingRate;
            //cout<<count<<endl;

            if(count>minsupport){

                Path temppath;
                temppath.pathidentifier=globalpathidcount;
                pathidentifier2pathsindex.insert({globalpathidcount,make_pair(length,paths[length].size())});
                globalpathidcount++;
                temppath.count=count;
                temppath.vertextuples=rulepath.first.vertextuples;
                temppath.edgelabels=rulepath.first.edgelabels;
                temppath.kleenestar=false;

                //temppath.Show();

                paths[length].push_back(temppath);

                #pragma omp parallel for
                for(int for_core=0;for_core<corenum;for_core++) {
                    for (auto matchedvertexid: matchedvertexsets[for_core]) {

                        matchedvertices[matchedvertexid].pathidentifiers[length].push_back(temppath.pathidentifier);
                        matchedvertices[matchedvertexid].targets[length].push_back(nexttargets[matchedvertexid]);
                    }
                }

                for(auto& pathidentifier: rulepath.second){
                    paths[length][pathidentifier2pathsindex[pathidentifier].second].horizontalExtendingPathidentifiers.push_back(temppath.pathidentifier);
                }


            }
        }
        //return;

        start=end;
        end=paths[length].size();
        numberOfCandidateAttribute++;
        if(start==end)break;
    }

}


void AssociationRules::ParallelFindFrequentLengthOnePath(OriginalGraph &graph, vector<Path>& paths0set, int maxattributesize){

    int test=0;

//    vector<vector<int>> parallelvertexidsets;
//    parallelvertexidsets.resize(corenum);
//    for(int i=0;i<graph.number_of_vertices;i++){
//        parallelvertexidsets[i%corenum].push_back(i);
//    }

    for(Path& path: paths0set){//extends paths, we call it path A

        //cout<<"test"<<endl;
        int count = 0;
        //vector<int> matchedvertexs;
        //matchedvertexs.clear();

        vector<vector<int>> matchedvertexids;
        matchedvertexids.resize(corenum);

//        matchedvertexids.clear();

        vector<vector<pair<int,int>>> matchedvertexsets;
        matchedvertexsets.resize(corenum);

        #pragma omp parallel for
        for(int for_core=0;for_core<corenum;for_core++){
            matchedvertexsets[for_core].clear();
        }

        vector<int> pathid(corenum);
        #pragma omp parallel for
        for(int for_core=0;for_core<corenum;for_core++) {
            matchedvertexsets[for_core]=ParallelPathContainment(parallelvertexidsets[for_core], 0, path.pathidentifier);
        }


        vector<bool> insertedtemptaregets;
        insertedtemptaregets.resize(graph.number_of_vertices);

        for(auto for_edgelabel: candidateedgelabels[0]) {//for all edge labels

            vector<vector<int>> temptargets;
            vector<vector<vector<int>>> targetslist;
            targetslist.resize(corenum);

            #pragma omp parallel for
            for(int for_core=0;for_core<corenum;for_core++) {
                targetslist[for_core]=ParallelFindTarget(graph, matchedvertexsets[for_core], for_edgelabel, 1);
            }

            vector<vector<int>> nexttargets;
            nexttargets.resize(graph.number_of_vertices);
            //vector<bool> matchedvertexresult;
            //matchedvertexresult.resize(graph.number_of_vertices);

            vector<bool> insertedmatchedvertexs;
            insertedmatchedvertexs.resize(graph.number_of_vertices);

            for(auto path0: targettuples_edges[0][for_edgelabel]) {//for all frequent vertex tuples

                if( path0.second.size() > maxattributesize)continue;
                test++;

//                cout<<path.vertextuples[0][0]<<"-"<<for_edgelabel<<"-"<<path0.second[0]<<endl;


                count=0;
                for(int i=0;i<graph.number_of_vertices;i++){
                    nexttargets[i].clear();
                    insertedmatchedvertexs[i]=false;
                }


                #pragma omp parallel for
                for(int for_core=0;for_core<corenum;for_core++) {
                     matchedvertexids[for_core].clear();
                     matchedvertexids[for_core]=ParallelVertexMatching(graph, matchedvertexsets[for_core], targetslist[for_core], nexttargets, path0.first);
                }

                for(int for_core=0;for_core<corenum;for_core++) {
                    count+=matchedvertexids[for_core].size();
                }

                //cout<<"test count:= " <<count<<endl;
                count*=1/ApproximateSamplingRate;
                //cout<<"test aproximate count:= "<< count<<endl;

                if (count > minsupport) {

                    Path temppath;
                    temppath.pathidentifier = globalpathidcount;
                    pathidentifier2pathsindex.insert({globalpathidcount,make_pair(1,paths[1].size())});
                    globalpathidcount++;
                    temppath.count = count;
                    temppath.vertextuples = path.vertextuples;
                    temppath.vertextuples.push_back(path0.second);
                    temppath.edgelabels.push_back(for_edgelabel);
                    temppath.kleenestar = false;

                    paths[1].push_back(temppath);

                    #pragma omp parallel for
                    for(int for_core=0;for_core<corenum;for_core++) {
                        for (auto matchedvertexid: matchedvertexids[for_core]) {

                            matchedvertices[matchedvertexid].pathidentifiers[1].push_back(temppath.pathidentifier);
                            matchedvertices[matchedvertexid].targets[1].push_back(nexttargets[matchedvertexid]);
                        }
                    }
                    paths[0][pathidentifier2pathsindex[path.pathidentifier].second].verticalExtendingPathidentifiers.push_back(temppath.pathidentifier);

                }
            }
        }
    }

}


void AssociationRules::ParallelFindFrequentSingleAttribute(OriginalGraph& graph){

    struct timespec temp1_startTime, temp1_endTime;
    double temp1time=0;

    vector<int> countedgelabels;
    countedgelabels.resize(graph.cardinality_of_edgelabels);
    vector<vector<int>> parallelcountedgelabels;
    parallelcountedgelabels.resize(corenum);

    vector<bool> tempcandidateedges(graph.cardinality_of_edgelabels,false);
    vector<int> countinverseedgelabels;
    countinverseedgelabels.resize(graph.cardinality_of_edgelabels);

    vector<vector<bool>> tempcandidateinverseedges;
    tempcandidateinverseedges.resize(maximumpathlength);
    for(int i=0;i<maximumpathlength;i++){
        tempcandidateinverseedges[i].resize(graph.cardinality_of_edgelabels);
        for(int j=0;j<maximumpathlength;j++){
            tempcandidateinverseedges[i][j]=false;
        }
    }

    vector<int> countindegrees_edgelabels;
    countindegrees_edgelabels.resize(graph.cardinality_of_edgelabels);

    vector<vector<int>> parallelcountinverseedgelabels;
    parallelcountinverseedgelabels.resize(corenum);

    for(int for_core=0;for_core<corenum;for_core++) {
        parallelcountedgelabels[for_core].resize(graph.cardinality_of_edgelabels);
        parallelcountinverseedgelabels[for_core].resize(graph.cardinality_of_edgelabels);
    }

    if(EmptysourceFlag){// path of empty attribute (i.e., all vertices)

        #pragma omp parallel for
        for (int i = 0; i < graph.cardinality_of_edgelabels; i++) {
            countedgelabels[i] = 0;
        }

        Path temppath;
        temppath.pathidentifier = globalpathidcount;
        pathidentifier2pathsindex.insert({globalpathidcount,make_pair(0,paths[0].size())});
        globalpathidcount++;
        temppath.count = graph.number_of_vertices;
        vector<int> rulevertextuple;
        rulevertextuple.clear();// i.e., empty
        temppath.vertextuples.push_back(rulevertextuple);
        temppath.kleenestar = false;

        paths[0].push_back(temppath);

        vector<int> temptarget;

        for (auto& vertex: graph.vertex_list) { // all vertices

            temptarget.clear();
            matchedvertices[vertex.vertexid].pathidentifiers[0].push_back(temppath.pathidentifier);
            //temptarget[for_core].clear();
            temptarget.push_back(vertex.vertexid);
            matchedvertices[vertex.vertexid].targets[0].push_back(temptarget);

            for (int for_edgelabel = 0; for_edgelabel < graph.cardinality_of_edgelabels; for_edgelabel++) {
                if (graph.vertex_connectingcount_edgelabels[vertex.vertexid][for_edgelabel] > 0)countedgelabels[for_edgelabel]++;
            }
        }

        matchedvertexs_attributeset.push_back(graph.connecting_vertex);

        #pragma omp parallel for
        for (int for_edgelabel = 0; for_edgelabel < graph.cardinality_of_edgelabels; for_edgelabel++) {
            //for (int for_core = 0; for_core < corenum; for_core++) countedgelabels[for_edgelabel]+=parallelcountedgelabels[for_core][for_edgelabel];
            if (countedgelabels[for_edgelabel] > minsupport)tempcandidateedges[for_edgelabel] = true;
        }

    }



    for(int for_attribute = 0; for_attribute < graph.cardinality_of_vertextuples; for_attribute++){

//        cout<<graph.numbers_of_attributes[for_attribute]<<endl;
        if(graph.numbers_of_attributes[for_attribute]<=minsupport)continue;

        #pragma omp parallel for
        for (int i = 0; i < graph.cardinality_of_edgelabels; i++) {
            countedgelabels[i] = 0;
        }

        Path temppath;
        temppath.pathidentifier = globalpathidcount;
        pathidentifier2pathsindex.insert({globalpathidcount,make_pair(0,paths[0].size())});
        globalpathidcount++;
        temppath.count = graph.numbers_of_attributes[for_attribute];
        vector<int> rulevertextuple;
        rulevertextuple.push_back(for_attribute);
        temppath.vertextuples.push_back(rulevertextuple);
        temppath.kleenestar = false;
        frequentattributes.push_back(for_attribute);

        paths[0].push_back(temppath);

        vector<int> temptarget;

        //for (int for_core = 0; for_core < corenum; for_core++) {
        for (auto& matchedvertexid: graph.attribute_connecting_vertex[for_attribute]) {

            temptarget.clear();
            matchedvertices[matchedvertexid].pathidentifiers[0].push_back(temppath.pathidentifier);
            //temptarget[for_core].clear();
            temptarget.push_back(matchedvertexid);
            matchedvertices[matchedvertexid].targets[0].push_back(temptarget);

            for (int for_edgelabel = 0; for_edgelabel < graph.cardinality_of_edgelabels; for_edgelabel++) {
                if (graph.vertex_connectingcount_edgelabels[matchedvertexid][for_edgelabel] > 0)countedgelabels[for_edgelabel]++;
            }
        }

        matchedvertexs_attributeset.push_back(graph.attribute_connecting_vertex[for_attribute]);

        #pragma omp parallel for
        for (int for_edgelabel = 0; for_edgelabel < graph.cardinality_of_edgelabels; for_edgelabel++) {
            //for (int for_core = 0; for_core < corenum; for_core++) countedgelabels[for_edgelabel]+=parallelcountedgelabels[for_core][for_edgelabel];
            if (countedgelabels[for_edgelabel] > minsupport)tempcandidateedges[for_edgelabel] = true;
        }
       // }

    }


    clock_gettime(CLOCK_REALTIME, &temp1_startTime);

    for(int for_targetattribute = 0; for_targetattribute < graph.cardinality_of_vertextuples; for_targetattribute++) {

        //cout<<"max : "<<graph.attribute_max_inverseconnectingcount_edgelabels[for_targetattribute]<<", "<<graph.maximumIndegree_edgelabel<<","<<minsupport / (graph.maximumIndegree_edgelabel^maximumpathlength)<<endl;
        if(graph.attribute_max_inverseconnectingcount_edgelabels[for_targetattribute] <= minsupport / (graph.maximumIndegree_edgelabel^maximumpathlength))continue;

        for (int for_edgelabel = 0; for_edgelabel < graph.cardinality_of_edgelabels; for_edgelabel++) {

            if (graph.attribute_inverseconnectingcount_edgelabels[for_targetattribute][for_edgelabel]==0 || graph.attribute_inverseconnectingcount_edgelabels[for_targetattribute][for_edgelabel] <= minsupport / (graph.attribute_inverseconnectingcount_edgelabels[for_targetattribute][for_edgelabel]*(graph.maximumIndegree_edgelabel^(maximumpathlength-1))))continue;
            //cout<<"each edge label: "<<graph.attribute_inverseconnectingcount_edgelabels[for_targetattribute][for_edgelabel]<<endl;

            for (int for_length = 0; for_length < maximumpathlength; for_length++) {

                double x;
                if(!ExtentionFlag||!EfficientFlag)x=-1;
                else if(for_length==0)x=1;
                //else x=ApproximateCandidateRate*(graph.attribute_inverseconnectingcount_edgelabels[for_targetattribute][for_edgelabel]*(graph.maximumIndegree_edgelabel^(for_length-1)));
                else x=pow(graph.maximumIndegree_edgelabel,(ApproximateCandidateRate*(for_length)));
                //cout<<x<<endl;sa


                if (graph.attribute_inverseconnectingcount_edgelabels[for_targetattribute][for_edgelabel]!=0&&graph.attribute_inverseconnectingcount_edgelabels[for_targetattribute][for_edgelabel] > minsupport / x){
                    //cout<<globaltargetcount<<","<<for_edgelabel<<","<<for_length<<","<<countinverseedgelabels[for_edgelabel]<<endl;

                    //cout<<"each edge label: "<<graph.attribute_inverseconnectingcount_edgelabels[for_targetattribute][for_edgelabel]<<endl;

                    vector<int> temptuples;
                    temptuples.push_back(for_targetattribute);
                    for (int for_length2 = for_length; for_length2 < maximumpathlength; for_length2++) {
                        targettuples_edges[for_length2][for_edgelabel].push_back(make_pair(globaltargetcount, temptuples));
                        tempcandidateinverseedges[for_length2][for_edgelabel] = true;
                    }

                    //#pragma omp parallel for
                    //for (int for_core = 0; for_core < corenum; for_core++) {
                    int testcount=0;
                    //for(int for_sourceattribute = 0; for_sourceattribute < graph.cardinality_of_vertextuples; for_sourceattribute++) {
                    for (auto &matchedvertexid: graph.attribute_connecting_vertex[for_targetattribute]) {
                        //for (auto &edgeid : graph.vertex_connecting_edge_list[matchedvertexid]) {
                         //   if (graph.edge_list[edgeid].edgelabel == for_edgelabel && Contains(graph.vertex_list[graph.edge_list[edgeid].dst].vertextuple, for_targetattribute)) {
                                matchedvertices[matchedvertexid].targetidentifiers.push_back(globaltargetcount);
                                //if(for_sourceattribute==0&&for_edgelabel==2)testcount++;
                         //       break;
                        //    }
                        //}
                    }
                    //}

                    //cout<<"testcount:: "<<testcount<<endl;

                    //}
                    globaltargetcount++;
                    break;
                }
            }
        }
    }

    //cout<<"size: "<<targettuples_edges[0][2]

    for(int for_edgelabel=0;for_edgelabel<graph.cardinality_of_edgelabels;for_edgelabel++) {
        //cout<<for_edgelabel<<":"<<tempcandidateedges[for_edgelabel]<<","<<tempcandidateinverseedges[for_edgelabel]<<endl;
        for (int for_length = 0; for_length < maximumpathlength;for_length++){
            if (tempcandidateedges[for_edgelabel] && tempcandidateinverseedges[for_length][for_edgelabel]){
                for (int for_length2 = for_length; for_length2 < maximumpathlength; for_length2++) {
                    candidateedgelabels[for_length2].push_back(for_edgelabel);
                }
                break;
            }
            if (for_length>0 && tempcandidateinverseedges[for_length][for_edgelabel]){
                for (int for_length2 = for_length; for_length2 < maximumpathlength; for_length2++) {
                    candidateedgelabels[for_length2].push_back(for_edgelabel);
                }
                break;
            }
        }
    }
    clock_gettime(CLOCK_REALTIME, &temp1_endTime);
    temp1time+=(temp1_endTime.tv_sec - temp1_startTime.tv_sec) + (temp1_endTime.tv_nsec - temp1_startTime.tv_nsec) * pow(10, -9);

    cout << "parallel single attribute temp1 = "<<temp1time<<endl;

}


void AssociationRules::ParallelEfficientCombinePaths(OriginalGraph graph){

    int combinecount=0;

    google::dense_hash_map<string, bool> findPatterns;
    findPatterns.set_empty_key("-1");

//    vector<vector<int>> parallelvertexidsets;
//    parallelvertexidsets.resize(corenum);
//    for(int i=0;i<graph.number_of_vertices;i++){
//        parallelvertexidsets[i%corenum].push_back(i);
//    }


    //Combining unit paths//

    vector<vector<int>> matchedvertexsets;
    int count_s;

    matchedvertexsets.resize(corenum);
    double startTime, endTime;
    struct timespec process_startTime, process_endTime;

    clock_gettime(CLOCK_REALTIME, &process_startTime);

    #pragma omp parallel shared(matchedvertexsets, count_s, findPatterns)
    {
        struct timespec  thread_startTime, thread_endTime;
        int for_core=omp_get_thread_num();
        double time=0;

        for (int for_length = minimumRulePathLength; for_length < maximumpathlength + 1; for_length++) {
            #pragma omp barrier
            CheckTime();
            //cout<<"length "<<for_length<<endl;
            int pathsize=paths[for_length].size();

            for (int for_path1 = 0; for_path1 < pathsize; for_path1++) {
                #pragma omp barrier
                //startTime = omp_get_wtime();
                Path path1;
                path1 = paths[for_length][for_path1];
                int pathid1;
                pathid1 = path1.pathidentifier;

                for (int i = 0; i < for_length + 1; i++) {
                    if (path1.vertextuples[i].size() > 1) {
                        continue;
                    }
                }


                //cout<<"path1-" << for_path1<<"/"<<paths[for_length].size()<<endl;
                for (int for_path2 = for_path1 + 1; for_path2 < pathsize; for_path2++) {
                    #pragma omp barrier
                    clock_gettime(CLOCK_REALTIME, &thread_startTime);

                    Path path2;

                    int pathid2;
                    bool continueFlag;
                    //cout<<"path2-" << for_path2<<endl;
                    path2 = paths[for_length][for_path2];
                    pathid2 = path2.pathidentifier;

                    //cout <<"test1: "<< for_path1 << "," << for_path2 << endl;
                    continueFlag = false;
                    for (int i = 0; i < for_length + 1; i++) {
                        if (path2.vertextuples[i].size() > 1) {
                            continueFlag = true;
                            break;
                        }
                    }
                    //cout << "test2-"<<for_core<<"id1, id2= "<<pathid1<<" , "<<pathid2 << endl;
//#pragma omp barrier
                    //cout<<"test_id-"<<for_core<<"id1, id2= "<<pathid1<<" , "<<pathid2<<endl;
                    if (continueFlag){
                        clock_gettime(CLOCK_REALTIME, &thread_endTime);
                        //endTime = omp_get_wtime();
                        time+= (thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9);
                        if((thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9)<0)cout<<"  whattttttttt!!!"<<endl;

                        continue;
                    }

//                if (path1.vertextuples[0].size() > 1 || path1.vertextuples[1].size() > 1)continue;
//                if (path2.vertextuples[0].size() > 1 || path2.vertextuples[1].size() > 1)continue;

                    if (NonDominanceOutputFlag&&CheckPathDominance(path1, path2))continueFlag = true;
//#pragma omp barrier
                    if (continueFlag){
                        clock_gettime(CLOCK_REALTIME, &thread_endTime);
                        //endTime = omp_get_wtime();
                        time+= (thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9);
                        if((thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9)<0)cout<<"  whattttttttt!!!"<<endl;

                        continue;
                    }

                    //cout<<"test3-"<<for_core<<endl;


                    count_s = 0;

                    combinecount++;

//#pragma omp parallel for
//                    for (int for_core = 0; for_core < corenum; for_core++) {

                    //#pragma omp barrier
                    //cout<<"test4-"<<for_core<<","<<matchedvertexsets.size()<<endl;
                    //matchedvertexsets[for_core].clear();

                    matchedvertexsets[for_core] = PatternVertexMatching_forEfficient(graph, parallelvertexidsets[for_core], for_length, for_length, pathid1, pathid2);

                    #pragma omp barrier

                    #pragma omp single
                    {
                        for (int for_core1 = 0; for_core1 < corenum; for_core1++) {
                            count_s += matchedvertexsets[for_core1].size();
                        }
                        count_s*=1/ApproximateSamplingRate;
                    }
                    //cout<<"count "<<count_s<<endl;


                    #pragma omp barrier
                    if (count_s > minsupport) {
                        Pattern temppattern;


                            temppattern.patternidentifier = patterns.size();
                            temppattern.count = count_s;
                            temppattern.duplicateFlag = false;

                            int temppathid1 = pathid1;
                            int temppathid2 = pathid2;

                            if (path1 > path2)swap(temppathid1, temppathid2);

                            temppattern.pathidentifiers.push_back(temppathid1);
                            temppattern.pathidentifiers.push_back(temppathid2);
                        #pragma omp single
                            {
                            patterns.push_back(temppattern);
                        }
//#pragma omp parallel for
//                        for (int for_core = 0; for_core < corenum; for_core++) {

                        int for_core=omp_get_thread_num();
                        for (auto matchedvertexid: matchedvertexsets[for_core]) {
                            matchedvertices[matchedvertexid].patternidentifiers.push_back(temppattern.patternidentifier);
                        }

                        //}
                        #pragma omp single
                        {
                            string temppair = to_string(pathid1) + "-" + to_string(pathid2);
                            if (path1 > path2) {
                                temppair = to_string(pathid2) + "-" + to_string(pathid1);
                            }
                            findPatterns[temppair] = true;
                            //cout<<"test6-"<<for_core<<endl;
                        }

//                 for (auto matchedvertexid: matchedvertexids_s) {
//                     matchedvertices[matchedvertexid].patternidentifiers.push_back(temppattern.patternidentifier);
//                     //matchedvertices[matchedvertexid].targets[for_length].push_back(nexttargets[matchedvertexid]);
//                 }

                    }

                    clock_gettime(CLOCK_REALTIME, &thread_endTime);
                    //endTime = omp_get_wtime();
                    time+= (thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9);
                    if((thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9)<0)cout<<"  whattttttttt!!!"<<endl;

                }
//                clock_gettime(CLOCK_REALTIME, &thread_endTime);
//                //endTime = omp_get_wtime();
//                time+= (thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9);
//                if((thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9)<0)cout<<"  whattttttttt!!!"<<endl;


            }

            //cout<<"length done" << for_length<<" ,core"<<omp_get_thread_num()<<endl;
        }
        #pragma omp critical
        cout<<"core ="<< for_core <<" : "<<time<< "[sec]"<<endl;
    }

    cout<<"unit path combination done: size ="<<patterns.size()<<endl;
    clock_gettime(CLOCK_REALTIME, &process_endTime);
    double time= (process_endTime.tv_sec - process_startTime.tv_sec) + (process_endTime.tv_nsec - process_startTime.tv_nsec) * pow(10, -9);
    cout<<" unit path combination time = "<<(process_endTime.tv_sec - process_startTime.tv_sec) + (process_endTime.tv_nsec - process_startTime.tv_nsec) * pow(10, -9)<<endl;


    int start=0;
    int end=patterns.size();

    int round=0;

    //std::ofstream fout("debuging", ios::app);
    clock_gettime(CLOCK_REALTIME, &process_startTime);

    #pragma omp parallel shared(matchedvertexsets, count_s, findPatterns)
    {
        struct timespec  thread_startTime, thread_endTime;
        int for_core = omp_get_thread_num();
        double time=0;
        while (1) {
            CheckTime();
            round++;
            for (int for_pattern = start; for_pattern < end; for_pattern++) {

                Pattern pathpattern = patterns[for_pattern];

                int pathid1 = pathpattern.pathidentifiers[0];
                int pathid2 = pathpattern.pathidentifiers[1];

                int oripathid1 = pathid1;
                int oripathid2 = pathid2;

                pair<int, int> pathindex1 = pathidentifier2pathsindex[pathid1];
                pair<int, int> pathindex2 = pathidentifier2pathsindex[pathid2];

                vector<int> matchedvertexids_s;
                //vector<int> matchedvertexids_st;
                matchedvertexids_s.clear();
                //matchedvertexids_st.clear();

                Path path1 = paths[pathindex1.first][pathindex1.second];
                Path path2 = paths[pathindex2.first][pathindex2.second];

//            if(for_pattern==23812){
//                cout<<"original path"<<endl;
//                path1.Show();
//                path2.Show();
//            }



//            vector<pair<int,int>> frequentextendingids1;
//            vector<pair<int,int>> frequentextendingids2;
//            frequentextendingids1.clear();
//            frequentextendingids2.clear();

                combinecount = 0;

                //cout<<"pathid1 ="<<pathid1<<"("<<pathindex1.first<<","<<pathindex1.second<<"), ex size="<<path1.verticalExtendingPathidentifiers.size()<<endl;
                //cout<<"pathid2 ="<<pathid2<<"("<<pathindex2.first<<","<<pathindex2.second<<"), ex size="<<path2.verticalExtendingPathidentifiers.size()<<endl;

                //#pragma omp barrier
                #pragma omp barrier
                for (auto &extendingidentifiers: path1.verticalExtendingPathidentifiers) {
                    //cout<<"start =" << start<<", vex1 id ="<<extendingidentifiers<<endl;
                    #pragma omp barrier
                    combinecount++;
                    pair<int, int> expathid = pathidentifier2pathsindex[extendingidentifiers];
                    Path expath = paths[expathid.first][expathid.second];


//                if(for_pattern==23812){
//                cout<<"path1 vertical"<<endl;
//                path1.Show();
//                path2.Show();
//                expath.Show();
//                }

                    if (NonDominanceOutputFlag&&CheckPathDominance(expath, path2)) {
                        //cout<<"Dominance"<<endl;
                        continue;
                    } else {
                        // cout<<"Non Dominance"<<endl;
                    }

                    //Duplicate Check: Check whether there is the same path in found patterns.
                    string temppair = to_string(extendingidentifiers) + "-" + to_string(pathid2);
                    if (expath > path2) {
                        temppair = to_string(pathid2) + "-" + to_string(extendingidentifiers);
                    }
                    if (findPatterns.find(temppair) != findPatterns.end())continue;
                    //#pragma omp barrier

//                    #pragma omp single
//                    {
//                        findPatterns[temppair] = true;
//                    }

                    count_s = 0;
                    combinecount++;

                    clock_gettime(CLOCK_REALTIME, &thread_startTime);
                    matchedvertexsets[for_core] = PatternVertexMatching_forEfficient(graph, parallelvertexidsets[for_core], pathindex1.first + 1, pathindex2.first, extendingidentifiers, pathid2);
                    clock_gettime(CLOCK_REALTIME, &thread_endTime);
                    time+= (thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9);


                    #pragma omp barrier

                    #pragma omp single
                    {
                        findPatterns[temppair] = true;
                        for (int for_core1 = 0; for_core1 < corenum; for_core1++) {
                            count_s += matchedvertexsets[for_core1].size();
                        }
                        count_s*=1/ApproximateSamplingRate;
                    }

                    #pragma omp barrier
                    if (count_s > minsupport) {
                        //frequentextendingids1.push_back(make_pair(extendingidentifiers,pathindex1.first+1));

                        Pattern temppattern;
                        temppattern.patternidentifier = patterns.size();
                        temppattern.count = count_s;


                        int tempextendingidentifier = extendingidentifiers;
                        int temppathid = pathid2;
                        if (expath > path2)swap(tempextendingidentifier, temppathid);


                        temppattern.pathidentifiers.push_back(tempextendingidentifier);
                        temppattern.pathidentifiers.push_back(temppathid);

                        #pragma omp single
                        {
                            patterns.push_back(temppattern);
                        }

                        for (auto matchedvertexid: matchedvertexsets[for_core]) {
                            matchedvertices[matchedvertexid].patternidentifiers.push_back(temppattern.patternidentifier);
                        }

                    }
                }

                #pragma omp barrier
                for (auto &extendingidentifiers: path2.verticalExtendingPathidentifiers) {
                    #pragma omp barrier

                    //cout<<"start =" << start<<", vex2 id ="<<extendingidentifiers<<endl;
                    combinecount++;

                    pair<int, int> expathid = pathidentifier2pathsindex[extendingidentifiers];
                    Path expath = paths[expathid.first][expathid.second];

//                cout<<"path2 vertical"<<endl;
//                path1.Show();
//                path2.Show();
//                expath.Show();

                    if (NonDominanceOutputFlag&&CheckPathDominance(path1, expath)) {
                        //cout<<"Dominance"<<endl;
                        continue;
                    } else {
                        //cout<<"Non Dominance"<<endl;
                    }

                    //Duplicate Check: Check whether there is the same path in found patterns.

                    string temppair = to_string(extendingidentifiers) + "-" + to_string(pathid1);
                    if (expath > path1) {
                        temppair = to_string(pathid1) + "-" + to_string(extendingidentifiers);
                    }
                    if (findPatterns.find(temppair) != findPatterns.end())continue;
                    //#pragma omp barrier

//                    #pragma omp single
//                    {
//                        findPatterns[temppair] = true;
//                    }


                    count_s = 0;

                    combinecount++;
                    clock_gettime(CLOCK_REALTIME, &thread_startTime);
                    matchedvertexsets[for_core] = PatternVertexMatching_forEfficient(graph, parallelvertexidsets[for_core], pathindex1.first, pathindex2.first + 1, pathid1, extendingidentifiers);
                    clock_gettime(CLOCK_REALTIME, &thread_endTime);
                    time+= (thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9);


                    #pragma omp barrier

                    #pragma omp single
                    {
                        findPatterns[temppair] = true;
                        for (int for_core1 = 0; for_core1 < corenum; for_core1++) {
                            count_s += matchedvertexsets[for_core1].size();
                        }
                        count_s*=1/ApproximateSamplingRate;
                    }

                    #pragma omp barrier
                    if (count_s > minsupport) {
                        //frequentextendingids2.push_back(make_pair(extendingidentifiers,pathindex2.first+1));

                        Pattern temppattern;
                        temppattern.patternidentifier = patterns.size();
                        temppattern.count = count_s;


                        int tempextendingidentifier = extendingidentifiers;
                        int temppathid = pathid1;
                        if (expath > path1)swap(tempextendingidentifier, temppathid);
                        //if(expath>path1)swap(extendingidentifiers, pathid1);


                        temppattern.pathidentifiers.push_back(temppathid);
                        temppattern.pathidentifiers.push_back(tempextendingidentifier);
                        #pragma omp single
                        {
                            patterns.push_back(temppattern);
                        }


                        for (auto matchedvertexid: matchedvertexsets[for_core]) {
                            matchedvertices[matchedvertexid].patternidentifiers.push_back(temppattern.patternidentifier);
                        }

                    }
                }
                #pragma omp barrier
                for (auto &extendingidentifiers: path1.horizontalExtendingPathidentifiers) {
                    #pragma omp barrier
                    //cout<<"start =" << start<<", hex1 id ="<<extendingidentifiers<<endl;
                    combinecount++;
                    //frequentextendingids1.push_back(make_pair(extendingidentifiers,pathindex1.first));

                    pair<int, int> expathid = pathidentifier2pathsindex[extendingidentifiers];
                    Path expath = paths[expathid.first][expathid.second];

//                if(for_pattern==23812) {
//                    cout<<"path1 horizontal"<<endl;
//                    path1.Show();
//                    path2.Show();
//                    expath.Show();
//                }


                    if (NonDominanceOutputFlag&&CheckPathDominance(expath, path2)) {
                        //cout<<"Dominance"<<endl;
                        continue;
                    } else {
                        //cout<<"Non Dominance"<<endl;
                    }
                    //Duplicate Check: Check whether there is the same path in found patterns.
                    string temppair = to_string(extendingidentifiers) + "-" + to_string(pathid2);
                    if (expath > path2) {
                        temppair = to_string(pathid2) + "-" + to_string(extendingidentifiers);
                    }
                    if (findPatterns.find(temppair) != findPatterns.end())continue;
                    //#pragma omp barrier

//                    #pragma omp single
//                    {
//                        findPatterns[temppair] = true;
//                    }

                    count_s = 0;
                    combinecount++;
                    clock_gettime(CLOCK_REALTIME, &thread_startTime);
                    matchedvertexsets[for_core] = PatternVertexMatching_forEfficient(graph, parallelvertexidsets[for_core], pathindex1.first, pathindex2.first, extendingidentifiers, pathid2);
                    clock_gettime(CLOCK_REALTIME, &thread_endTime);
                    time+= (thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9);

                    #pragma omp barrier

                    #pragma omp single
                    {
                        findPatterns[temppair] = true;
                        for (int for_core1 = 0; for_core1 < corenum; for_core1++) {
                            count_s += matchedvertexsets[for_core1].size();
                        }
                        count_s*=1/ApproximateSamplingRate;
                    }

                    #pragma omp barrier
                    if (count_s > minsupport) {

                        //frequentextendingids1.push_back(make_pair(extendingidentifiers,pathindex1.first));

                        Pattern temppattern;
                        temppattern.patternidentifier = patterns.size();
                        temppattern.count = count_s;


                        int tempextendingidentifier = extendingidentifiers;
                        int temppathid = pathid2;
                        if (expath > path2)swap(tempextendingidentifier, temppathid);
                        //if(expath>path2)swap(extendingidentifiers, pathid2);


                        temppattern.pathidentifiers.push_back(tempextendingidentifier);
                        temppattern.pathidentifiers.push_back(temppathid);


                        #pragma omp single
                        {
                            patterns.push_back(temppattern);
                        }

                        for (auto matchedvertexid: matchedvertexsets[for_core]) {
                            matchedvertices[matchedvertexid].patternidentifiers.push_back(temppattern.patternidentifier);
                        }
                    }
                }

                #pragma omp barrier
                for (auto &extendingidentifiers: path2.horizontalExtendingPathidentifiers) {
                    #pragma omp barrier

                    //cout<<"start =" << start<<", hex2 id ="<<extendingidentifiers<<endl;
                    combinecount++;
                    //frequentextendingids2.push_back(make_pair(extendingidentifiers,pathindex2.first));
                    pair<int, int> expathid = pathidentifier2pathsindex[extendingidentifiers];
                    Path expath = paths[expathid.first][expathid.second];

//                cout<<"path2 horizontal"<<endl;
//                path1.Show();
//                path2.Show();
//                expath.Show();
                    startTime = omp_get_wtime();
                    if (NonDominanceOutputFlag&&CheckPathDominance(path1, expath)) {
                        //cout<<"Dominance"<<endl;
                        continue;
                    } else {
                        //cout<<"Non Dominance"<<endl;
                    }

                    //Duplicate Check: Check whether there is the same path in found patterns.
                    string temppair = to_string(extendingidentifiers) + "-" + to_string(pathid1);
                    if (expath > path1) {
                        temppair = to_string(pathid1) + "-" + to_string(extendingidentifiers);
                    }
                    if (findPatterns.find(temppair) != findPatterns.end())continue;
                    //#pragma omp barrier
//                    #pragma omp single
//                    {
//                        findPatterns[temppair] = true;
//                    }

                    count_s = 0;
                    combinecount++;

                    clock_gettime(CLOCK_REALTIME, &thread_startTime);
                    matchedvertexsets[for_core] = PatternVertexMatching_forEfficient(graph, parallelvertexidsets[for_core], pathindex1.first, pathindex2.first, pathid1, extendingidentifiers);
                    clock_gettime(CLOCK_REALTIME, &thread_endTime);
                    time+= (thread_endTime.tv_sec - thread_startTime.tv_sec) + (thread_endTime.tv_nsec - thread_startTime.tv_nsec) * pow(10, -9);

                    #pragma omp barrier

                    #pragma omp single
                    {
                        findPatterns[temppair] = true;
                        for (int for_core1 = 0; for_core1 < corenum; for_core1++) {
                            count_s += matchedvertexsets[for_core1].size();
                        }
                        count_s*=1/ApproximateSamplingRate;
                    }


                    #pragma omp barrier
                    if (count_s > minsupport) {

                        //frequentextendingids2.push_back(make_pair(extendingidentifiers,pathindex2.first));

                        Pattern temppattern;
                        temppattern.patternidentifier = patterns.size();
                        temppattern.count = count_s;

                        int tempextendingidentifier = extendingidentifiers;
                        int temppathid = pathid1;
                        if (expath > path1)swap(tempextendingidentifier, temppathid);

                        temppattern.pathidentifiers.push_back(temppathid);
                        temppattern.pathidentifiers.push_back(tempextendingidentifier);

                        #pragma omp single
                        {
                            patterns.push_back(temppattern);
                        }

                        for (auto matchedvertexid: matchedvertexsets[for_core]) {
                            matchedvertices[matchedvertexid].patternidentifiers.push_back(temppattern.patternidentifier);
                        }
                    }
                }
            }

            start = end;
            end = patterns.size();
            if (start == end)break;
            //if(round==2)break;
            cout << "total size = " << patterns.size() << endl;
        }

        #pragma omp critical
        cout<<"core ="<< for_core <<" : "<<time<< "[sec]"<<endl;

        //cout<<combinecount<<endl;
    }

    clock_gettime(CLOCK_REALTIME, &process_endTime);
    cout<<" combine path find time = "<<(process_endTime.tv_sec - process_startTime.tv_sec) + (process_endTime.tv_nsec - process_startTime.tv_nsec) * pow(10, -9)<<endl;


}


void AssociationRules::ParallelCombinePaths(OriginalGraph graph){

//    vector<vector<int>> parallelvertexidsets;
//    parallelvertexidsets.resize(corenum);
//    for(int i=0;i<graph.number_of_vertices;i++){
//        parallelvertexidsets[i%corenum].push_back(i);
//    }

    for(int for_length1=minimumRulePathLength;for_length1<maximumpathlength+1;for_length1++){

        for(int for_length2=for_length1;for_length2<maximumpathlength+1;for_length2++){
            CheckTime();
            for(int for_path1 = 0; for_path1 < paths[for_length1].size();for_path1++){

                int start =0;
                if(for_length1==for_length2) start = for_path1+1;
                else start=0;

                for(int for_path2 = start; for_path2 < paths[for_length2].size();for_path2++){

                     //cout<<for_path1<<","<<for_path2<<endl;

                     Path path1 =  paths[for_length1][for_path1];
                     Path path2 =  paths[for_length2][for_path2];
                     int pathid1 = path1.pathidentifier;
                     int pathid2 = path2.pathidentifier;

                     vector<int> matchedvertexids_s;
                     vector<int> matchedvertexids_st;
                     matchedvertexids_s.clear();
                     //matchedvertexids_st.clear();

                     //path1.Show();
                     //path2.Show();

                     if (NonDominanceOutputFlag&&CheckPathDominance(path1, path2)){
                         //cout<<"Dominance"<<endl;
                         continue;
                     }
                     //else cout<< "Non Dominance"<<endl;

                     vector<vector<int>> matchedvertexsets;
                     matchedvertexsets.resize(corenum);
                     int count_s = 0;

                     #pragma omp parallel for
                     for (int for_core = 0; for_core < corenum; for_core++) {
                         matchedvertexsets[for_core] = PatternVertexMatching_forEfficient(graph, parallelvertexidsets[for_core], for_length1, for_length2, pathid1, pathid2);
                     }

                     for (int for_core = 0; for_core < corenum; for_core++) {
                        count_s += matchedvertexsets[for_core].size();
                     }


                     if(count_s>minsupport){
                         Pattern temppattern;
                         temppattern.patternidentifier = patterns.size();
                         temppattern.count = count_s;

                         if(path2>path1)swap(pathid1, pathid2);

                         temppattern.pathidentifiers.push_back(pathid1);
                         temppattern.pathidentifiers.push_back(pathid2);

                         patterns.push_back(temppattern);
                         #pragma omp parallel for
                         for (int for_core = 0; for_core < corenum; for_core++) {
                            for (auto matchedvertexid: matchedvertexsets[for_core]) {
                                matchedvertices[matchedvertexid].patternidentifiers.push_back(temppattern.patternidentifier);
                            }
                         }

                     }

                }
            }
        }
    }
}




//////////////////////////////////////////////////////



bool AssociationRules::VertexAttributeMatching(vector<int> &graphvertexattributes_, vector<int> &rulevertexattributes_){
//Efficientcy can be improved


    for (int for_attribute = 0; for_attribute < rulevertexattributes_.size(); for_attribute++) {

//        if (rulevertexattributes_[for_attribute] == -1)continue;
        if (!Contains(graphvertexattributes_,rulevertexattributes_[for_attribute])) {
            return false;
        }
    }

    return true;

//    if(notzeroFlag)return true;
//    else return false;
}


bool AssociationRules::CheckVertexTupleDominance(vector<int>& vertextuple1, vector<int>& vertextuple2){
//if vertex1 dominates vertex2, return true;


    if(vertextuple1.size() > vertextuple2.size())return false;

    int tuple1pos=0;
    int tuple2pos=0;
    while(1){

        if(tuple1pos==vertextuple1.size() || tuple2pos==vertextuple2.size())break;

        if (vertextuple1[tuple1pos] == vertextuple2[tuple2pos]){
            tuple1pos++;
            tuple2pos++;
        }
        else if(vertextuple1[tuple1pos]<vertextuple2[tuple2pos]) {
            tuple1pos++;
            if(tuple1pos==vertextuple1.size())return false;
        }
        else return false;

    }

    //if(tuple2pos!=vertextuple2.size())return false;

    return true;
}

bool AssociationRules::CheckPathDominance(Path& path1, Path& path2){
//if path1 (resp. path2) dominates path2 (resp. path1), return true;

    int length=0;
    bool dominancepath1=false;
    bool dominancepath2=false;

    // pick longer length of paths
    if(path1.vertextuples.size() < path2.vertextuples.size()){
        length = path1.vertextuples.size();
    }
    else {
        length = path2.vertextuples.size();
    }

    //check edge labe deffierence. If one of the edge labels are different, they are not dominated each other//
    for(int for_length=0;for_length<length; for_length++) {
        if (for_length < length - 1) {
            if (path1.edgelabels[for_length] != path2.edgelabels[for_length]) {
                return false;
            }
        }
    }

    for(int for_length=0;for_length<length; for_length++){

        int tuple1pos=0;
        int tuple2pos=0;

        if(path1.vertextuples[for_length].size()==path2.vertextuples[for_length].size()){
            for(int for_tuple=0;for_tuple<path1.vertextuples[for_length].size();for_tuple++)
            {
                if(path1.vertextuples[for_length][for_tuple] != path2.vertextuples[for_length][for_tuple])return false;
            }
        }
        else {

            while (1) {

                if (tuple1pos == path1.vertextuples[for_length].size() && tuple2pos == path2.vertextuples[for_length].size())break;
                else if (tuple1pos == path1.vertextuples[for_length].size()) {
                    dominancepath1 = true;
                    if (dominancepath2)return false;
                    else break;
                } else if (tuple2pos == path2.vertextuples[for_length].size()) {
                    dominancepath2 = true;
                    if (dominancepath1)return false;
                    else break;
                }

                if (path1.vertextuples[for_length][tuple1pos] == path2.vertextuples[for_length][tuple2pos]) {
                    tuple1pos++;
                    tuple2pos++;
                } else if (path1.vertextuples[for_length][tuple1pos] < path2.vertextuples[for_length][tuple2pos]) {
                    tuple1pos++;
                    dominancepath2 = true;
                    if (dominancepath1)return false;
                } else {
                    tuple2pos++;
                    dominancepath1 = true;
                    if (dominancepath2)return false;

                }
            }
        }

        if(dominancepath1&&path1.vertextuples.size() > path2.vertextuples.size())return false;
        if(dominancepath2&&path2.vertextuples.size() > path1.vertextuples.size())return false;
    }

    return true;
}




Path AssociationRules::MergeTwoPath(Path& path1, Path& path2, int numberOfCandidateAttribute) {

    Path outputpath;
    outputpath.clear();

    if(path1.kleenestar!=path2.kleenestar)return outputpath;
    if(path1.edgelabels.size()!=path2.edgelabels.size())return outputpath;
    else{
        for(int for_edgelabel=0;for_edgelabel<path1.edgelabels.size();for_edgelabel++){
            if(path1.edgelabels[for_edgelabel]!=path2.edgelabels[for_edgelabel])return outputpath;
            outputpath.edgelabels.push_back(path1.edgelabels[for_edgelabel]);
        }
    }

    for(int for_tuples=0;for_tuples<path1.vertextuples.size();for_tuples++){

        vector<int> tempvertextuple;
        tempvertextuple.clear();
        int path1_att = 0;
        int path2_att = 0;

        int numberOfMathcedAttribute = 0;
        bool diffFlag = true;

        for (int for_att = 0; for_att < numberOfCandidateAttribute - 1; for_att++) {
            diffFlag = false;
            if (path1.vertextuples[for_tuples][for_att] != path2.vertextuples[for_tuples][for_att]) {
                numberOfCandidateAttribute = for_att;
                diffFlag = true;
                break;
            }
        }

        if (numberOfMathcedAttribute == numberOfCandidateAttribute - 1 && diffFlag) {
            vector<int> candidatetuple;
            for (int for_att = 0; for_att < numberOfCandidateAttribute - 1; for_att++) {
                candidatetuple.push_back(path1.vertextuples[for_tuples][for_att]);
            }

            int att1 = path1.vertextuples[for_tuples][numberOfCandidateAttribute - 1];
            int att2 = path2.vertextuples[for_tuples][numberOfCandidateAttribute - 1];

            if (att1 == att2) {
                outputpath.clear();
                return outputpath;
            }
            if (att1 > att2)swap(att1, att2);

            candidatetuple.push_back(att1), candidatetuple.push_back(att2);
            outputpath.vertextuples.push_back(candidatetuple);
        } else {
            outputpath.clear();
            return outputpath;
        }

    }

    outputpath.kleenestar=path1.kleenestar;
    return outputpath;

}


Path AssociationRules::MergeTwoPathGivenTuples(Path& path1, Path& path2, int givenNumberOfCandidateAttribute, vector<int>& mergetuple) {

    Path outputpath;
    outputpath.clear();

    if(path1.kleenestar!=path2.kleenestar)return outputpath;
    if(path1.edgelabels.size()!=path2.edgelabels.size())return outputpath;
    else{
        for(int for_edgelabel=0;for_edgelabel<path1.edgelabels.size();for_edgelabel++){
            if(path1.edgelabels[for_edgelabel]!=path2.edgelabels[for_edgelabel])return outputpath;
            outputpath.edgelabels.push_back(path1.edgelabels[for_edgelabel]);
        }
    }

    for(int for_tuples=0;for_tuples<path1.vertextuples.size();for_tuples++){

        vector<int> tempvertextuple;
        tempvertextuple.clear();
        int path1_att = 0;
        int path2_att = 0;

        int numberOfMathcedAttribute = 0;
        bool diffFlag = true;

        if(!path1.vertextuples[for_tuples].empty()&&!path2.vertextuples[for_tuples].empty()&&Contains(mergetuple, for_tuples)){

            //cout<<"test1 MergeTwoPathGivenTuples: "<<for_tuples<<endl;

            int numberOfCandidateAttribute=givenNumberOfCandidateAttribute;

            for (int for_att = 0; for_att < numberOfCandidateAttribute - 1; for_att++) {
                diffFlag = false;
                if (path1.vertextuples[for_tuples][for_att] != path2.vertextuples[for_tuples][for_att]) {
                    numberOfCandidateAttribute = for_att;
                    diffFlag = true;
                    break;
                }
            }

            if (numberOfMathcedAttribute == numberOfCandidateAttribute - 1 && diffFlag) {
                vector<int> candidatetuple;
                for (int for_att = 0; for_att < numberOfCandidateAttribute - 1; for_att++) {
                    candidatetuple.push_back(path1.vertextuples[for_tuples][for_att]);
                }

                int att1 = path1.vertextuples[for_tuples][numberOfCandidateAttribute - 1];
                int att2 = path2.vertextuples[for_tuples][numberOfCandidateAttribute - 1];

                if (att1 == att2) {
                    outputpath.clear();
                    return outputpath;
                }
                if (att1 > att2)swap(att1, att2);

                candidatetuple.push_back(att1), candidatetuple.push_back(att2);
                outputpath.vertextuples.push_back(candidatetuple);
            } else {
                outputpath.clear();
                return outputpath;
            }
        }
        else{
            //cout<<"test2 MergeTwoPathGivenTuples: "<<for_tuples<<endl;

            if(path1.vertextuples[for_tuples].size()!=path2.vertextuples[for_tuples].size()){
                outputpath.clear();
                return outputpath;
            }
            for (int for_att = 0; for_att < path1.vertextuples[for_tuples].size(); for_att++) {
                if (path1.vertextuples[for_tuples][for_att] != path2.vertextuples[for_tuples][for_att]) {
                    outputpath.clear();
                    return outputpath;
                }
            }

            outputpath.vertextuples.push_back(path1.vertextuples[for_tuples]);

        }

    }

    outputpath.kleenestar=path1.kleenestar;
    return outputpath;

}




Path AssociationRules::MergeTwoPathSourceAndTarget(Path& path1, Path& path2, int numberOfCandidateAttribute) {

    Path outputpath;
    outputpath.clear();

    if(path1.kleenestar!=path2.kleenestar)return outputpath;
    if(path1.edgelabels.size()!=path2.edgelabels.size())return outputpath;
    else{
        for(int for_edgelabel=0;for_edgelabel<path1.edgelabels.size();for_edgelabel++){
            if(path1.edgelabels[for_edgelabel]!=path2.edgelabels[for_edgelabel])return outputpath;
            outputpath.edgelabels.push_back(path1.edgelabels[for_edgelabel]);
        }
    }

    for(int for_tuples=0;for_tuples<path1.vertextuples.size();for_tuples++){

        if(for_tuples==0 || for_tuples==path1.vertextuples.size()-1) {
            vector<int> tempvertextuple;
            tempvertextuple.clear();
            int path1_att = 0;
            int path2_att = 0;

            int numberOfMathcedAttribute = 0;
            bool diffFlag = true;

            for (int for_att = 0; for_att < numberOfCandidateAttribute - 1; for_att++) {
                diffFlag = false;
                if (path1.vertextuples[for_tuples][for_att] != path2.vertextuples[for_tuples][for_att]) {
                    numberOfCandidateAttribute = for_att;
                    diffFlag = true;
                    break;
                }
            }

            if (numberOfMathcedAttribute == numberOfCandidateAttribute - 1 && diffFlag) {
                vector<int> candidatetuple;
                for (int for_att = 0; for_att < numberOfCandidateAttribute - 1; for_att++) {
                    candidatetuple.push_back(path1.vertextuples[for_tuples][for_att]);
                }

                int att1 = path1.vertextuples[for_tuples][numberOfCandidateAttribute - 1];
                int att2 = path2.vertextuples[for_tuples][numberOfCandidateAttribute - 1];

                if (att1 == att2) {
                    outputpath.clear();
                    return outputpath;
                }
                if (att1 > att2)swap(att1, att2);

                candidatetuple.push_back(att1), candidatetuple.push_back(att2);
                outputpath.vertextuples.push_back(candidatetuple);
            } else {
                outputpath.clear();
                return outputpath;
            }
        }
        else{
            outputpath.vertextuples.push_back(path1.vertextuples[for_tuples]);
        }

    }

    outputpath.kleenestar=path1.kleenestar;
    return outputpath;

}



Path AssociationRules::MergeTwoOneLengthPathSource(Path& path1, Path& path2, int numberOfCandidateAttribute) {

    Path outputpath;
    outputpath.clear();


    if(path1.kleenestar!=path2.kleenestar)return outputpath;

    if(path1.edgelabels.size()!=path2.edgelabels.size())return outputpath;
    else{
        for(int for_edgelabel=0;for_edgelabel<path1.edgelabels.size();for_edgelabel++){
            if(path1.edgelabels[for_edgelabel]!=path2.edgelabels[for_edgelabel])return outputpath;
            outputpath.edgelabels.push_back(path1.edgelabels[for_edgelabel]);
        }
    }


    vector<int> tempvertextuple;
    tempvertextuple.clear();
    int path1_att=0;
    int path2_att=0;

    int numberOfMathcedAttribute=0;
    bool diffFlag=true;
    for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
        diffFlag=false;
        if(path1.vertextuples[0][for_att] != path2.vertextuples[0][for_att]) {
            numberOfCandidateAttribute=for_att;
            diffFlag=true;
            break;
        }
    }

    if(numberOfMathcedAttribute == numberOfCandidateAttribute-1&&diffFlag){
        vector<int> candidatetuple;
        for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
            candidatetuple.push_back(path1.vertextuples[0][for_att]);
        }

        int att1 = path1.vertextuples[0][numberOfCandidateAttribute-1];
        int att2 = path2.vertextuples[0][numberOfCandidateAttribute-1];
        if(att1==att2){
            outputpath.clear();
            return outputpath;
        }
        if(att1>att2)swap(att1,att2);

        candidatetuple.push_back(att1),candidatetuple.push_back(att2);
        outputpath.vertextuples.push_back(candidatetuple);
    }
    else{
        outputpath.clear();
        return outputpath;
    }


    outputpath.vertextuples.push_back(path1.vertextuples[1]);
    outputpath.kleenestar=path1.kleenestar;
    return outputpath;

}



Path AssociationRules::MergeTwoOneLengthPathTarget(Path& path1, Path& path2, int numberOfCandidateAttribute) {

    Path outputpath;
    outputpath.clear();

    if(path1.kleenestar!=path2.kleenestar)return outputpath;

    if(path1.edgelabels.size()!=path2.edgelabels.size())return outputpath;
    else{
        for(int for_edgelabel=0;for_edgelabel<path1.edgelabels.size();for_edgelabel++){
            if(path1.edgelabels[for_edgelabel]!=path2.edgelabels[for_edgelabel])return outputpath;
            outputpath.edgelabels.push_back(path1.edgelabels[for_edgelabel]);
        }
    }

    outputpath.vertextuples.push_back(path1.vertextuples[0]);

    vector<int> tempvertextuple;
    tempvertextuple.clear();
    int path1_att=0;
    int path2_att=0;

    int numberOfMathcedAttribute=0;
    bool diffFlag=true;
    for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
        diffFlag=false;
        if(path1.vertextuples[1][for_att] != path2.vertextuples[1][for_att]) {
            numberOfCandidateAttribute=for_att;
            diffFlag=true;
            break;
        }
    }

    if(numberOfMathcedAttribute == numberOfCandidateAttribute-1&&diffFlag){
        vector<int> candidatetuple;
        for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
            candidatetuple.push_back(path1.vertextuples[1][for_att]);
        }

        int att1 = path1.vertextuples[1][numberOfCandidateAttribute-1];
        int att2 = path2.vertextuples[1][numberOfCandidateAttribute-1];
        if(att1==att2){
            outputpath.clear();
            return outputpath;
        }
        if(att1>att2)swap(att1,att2);

        candidatetuple.push_back(att1),candidatetuple.push_back(att2);
         outputpath.vertextuples.push_back(candidatetuple);
    }
    else{
        outputpath.clear();
        return outputpath;
    }

    outputpath.kleenestar=path1.kleenestar;
    return outputpath;

}


Path AssociationRules::MergeTwoLLengthPathSource(Path& path1, Path& path2, int length, int numberOfCandidateAttribute) {

    Path outputpath;
    outputpath.clear();

    //Edge check//
    if(path1.edgelabels.size()!=path2.edgelabels.size())return outputpath;
    else{
        for(int for_edgelabel=0;for_edgelabel<path1.edgelabels.size();for_edgelabel++){
            if(path1.edgelabels[for_edgelabel]!=path2.edgelabels[for_edgelabel])return outputpath;
            outputpath.edgelabels.push_back(path1.edgelabels[for_edgelabel]);
        }
    }


    for(int for_length=1; for_length < length+1; for_length++){

        if(path1.vertextuples[for_length].size()!=path2.vertextuples[for_length].size())return outputpath;

        for(int for_att=0;for_att<path1.vertextuples[for_length].size(); for_att++) {
            if(path1.vertextuples[for_length][for_att] != path2.vertextuples[for_length][for_att]) {
                return outputpath;
            }
        }
    }


    vector<int> tempvertextuple;
    tempvertextuple.clear();
    int path1_att=0;
    int path2_att=0;


    int numberOfMathcedAttribute=0;
    bool diffFlag=true;
    for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
        diffFlag=false;
        if(path1.vertextuples[0][for_att] != path2.vertextuples[0][for_att]) {
            numberOfCandidateAttribute=for_att;
            diffFlag=true;
            break;
        }
    }

    if(numberOfMathcedAttribute == numberOfCandidateAttribute-1&&diffFlag){
        vector<int> candidatetuple;
        for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
            candidatetuple.push_back(path1.vertextuples[0][for_att]);
        }

        int att1 = path1.vertextuples[0][numberOfCandidateAttribute-1];
        int att2 = path2.vertextuples[0][numberOfCandidateAttribute-1];
        if(att1==att2){
            outputpath.clear();
            return outputpath;
        }
        if(att1>att2)swap(att1,att2);

        candidatetuple.push_back(att1),candidatetuple.push_back(att2);
         outputpath.vertextuples.push_back(candidatetuple);
    }
    else{
        outputpath.clear();
        return outputpath;
    }

    for(int for_length=1; for_length < length+1; for_length++) {
        outputpath.vertextuples.push_back(path1.vertextuples[for_length]);
    }

    outputpath.kleenestar=false;
    return outputpath;

}


Path AssociationRules::MergeTwoLLengthPathTarget(Path& path1, Path& path2, int length, int numberOfCandidateAttribute) {

    Path outputpath;
    outputpath.clear();

    if(path1.edgelabels.size()!=path2.edgelabels.size())return outputpath;
    else{
        for(int for_edgelabel=0;for_edgelabel<path1.edgelabels.size();for_edgelabel++){
            if(path1.edgelabels[for_edgelabel]!=path2.edgelabels[for_edgelabel])return outputpath;
            outputpath.edgelabels.push_back(path1.edgelabels[for_edgelabel]);
        }
    }

    for(int for_length=0; for_length < length; for_length++){

        if(path1.vertextuples[for_length].size()!=path2.vertextuples[for_length].size())return outputpath;

        for(int for_att=0;for_att<path1.vertextuples[for_length].size(); for_att++) {
            if(path1.vertextuples[for_length][for_att] != path2.vertextuples[for_length][for_att]) {
                return outputpath;
            }
        }
    }

    for(int for_length=0; for_length < length; for_length++) {
        outputpath.vertextuples.push_back(path1.vertextuples[for_length]);
    }
    vector<int> tempvertextuple;
    tempvertextuple.clear();
    int path1_att=0;
    int path2_att=0;


    int numberOfMathcedAttribute=0;
    bool diffFlag=true;
    for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
        diffFlag=false;
        if(path1.vertextuples[length][for_att] != path2.vertextuples[length][for_att]) {
            numberOfCandidateAttribute=for_att;
            diffFlag=true;
            break;
        }
    }

    if(numberOfMathcedAttribute == numberOfCandidateAttribute-1&&diffFlag){
        vector<int> candidatetuple;
        for(int for_att=0;for_att<numberOfCandidateAttribute-1; for_att++) {
            candidatetuple.push_back(path1.vertextuples[length][for_att]);
        }

        int att1 = path1.vertextuples[length][numberOfCandidateAttribute-1];
        int att2 = path2.vertextuples[length][numberOfCandidateAttribute-1];
        if(att1==att2){
            outputpath.clear();
            return outputpath;
        }
        if(att1>att2)swap(att1,att2);

        candidatetuple.push_back(att1),candidatetuple.push_back(att2);
         outputpath.vertextuples.push_back(candidatetuple);
    }
    else{
        outputpath.clear();
        return outputpath;
    }

    outputpath.kleenestar=false;
    return outputpath;

}





void AssociationRules::ParallelVertexAttributeMatching(OriginalGraph& graph, vector<int>& matchedvertexids, vector<int>& vertexidset, int attribute, vector<int> & countinverseedgelabels){

    //vector<int> matchedvertexids;
    //matchedvertexids.clear();

    for(auto vertexid: vertexidset){

        if(Contains(graph.vertex_list[vertexid].vertextuple,attribute)){
            matchedvertexids.push_back(vertexid);

            for(int for_edgelabel=0;for_edgelabel<graph.cardinality_of_edgelabels;for_edgelabel++){
                countinverseedgelabels[for_edgelabel]+=graph.vertex_inverseconnectingcount_edgelabels[vertexid][for_edgelabel];
            }

        }
    }

    //return matchedvertexids;
}



vector<int> AssociationRules::ParallelVertexMatching(OriginalGraph & graph, vector<pair<int,int>>& vertexidset, vector<vector<int>>& targetslist, vector<vector<int>>& nexttargets, int path0id){

    vector<int> matchedvertexids;
    matchedvertexids.clear();
    vector<bool> insertedmatchedvertexs;
    insertedmatchedvertexs.resize(graph.number_of_vertices);

    for (int for_vertex=0; for_vertex < vertexidset.size(); for_vertex++) {
        //cout<<for_vertex<<","<<numberOfmatchedvertex<<endl;
        for (auto target: targetslist[for_vertex]) {

            if (Contains(matchedvertices[target].targetidentifiers,path0id)) {

                if (!insertedmatchedvertexs[vertexidset[for_vertex].first]) {
                    matchedvertexids.push_back(vertexidset[for_vertex].first);
                    insertedmatchedvertexs[vertexidset[for_vertex].first] = true;
                }
                nexttargets[vertexidset[for_vertex].first].push_back(target);
            }
        }
    }
    return matchedvertexids;

}


vector<int> AssociationRules::ParallelVertexMultipleAttributeMatching(OriginalGraph& graph, vector<int>& vertexidset, vector<int>& pathidset){

    vector<int> matchedvertexids;
    int count=0;

    for (auto vertexid: vertexidset) {

        bool matchFlag = true;
        bool firstFlag = true;

        for (int for_attribute = 0; for_attribute < pathidset.size(); for_attribute++) {

//        if (rulevertexattributes_[for_attribute] == -1)continue;
            if (!Contains(graph.vertex_list[vertexid].vertextuple, pathidset[for_attribute])) {
                matchFlag = false;
                break;
            }
        }
        if (matchFlag) {
            matchedvertexids.push_back(vertexid);
        }
    }

    return matchedvertexids;

}



vector<int> AssociationRules::ParallelVertexMatching_forEfficient(OriginalGraph & graph, vector<int>& vertexidset, vector<int>& pathidset, int length, vector<vector<int>>& nexttargets){

    vector<int> matchedvertexids;
    int count=0;

    for (auto vertexid: vertexidset) {

        bool matchFlag=true;
        bool firstFlag=true;

        for(auto& pathidentifier: pathidset){
            int pathid = Contains_ReturnIndex(matchedvertices[vertexid].pathidentifiers[length], pathidentifier);//vertex matching the path A
            if (pathid == -1){
                matchFlag=false;
                break;
            }
            if(firstFlag){
                nexttargets[vertexid]=matchedvertices[vertexid].targets[length][pathid];
                firstFlag=false;
            }
            else{
                nexttargets[vertexid] = intersection(nexttargets[vertexid], matchedvertices[vertexid].targets[length][pathid]);

                if(nexttargets[vertexid].empty()){
                    matchFlag=false;
                    break;
                }
            }
        }

        if(matchFlag){
            matchedvertexids.push_back(vertexid);
        }
    }

    return matchedvertexids;

}


vector<int> AssociationRules::PatternVertexMatching_forEfficient(OriginalGraph &, vector<int>& vertexidset, int length1, int length2, int pathid1, int pathid2){

    vector<int> matchedvertexids;
    for (auto vertexid: vertexidset) {

        int pathid1_index = Contains_ReturnIndex(matchedvertices[vertexid].pathidentifiers[length1], pathid1);
        int pathid2_index = Contains_ReturnIndex(matchedvertices[vertexid].pathidentifiers[length2], pathid2);

//        if(pathid1_index>-1 || pathid2_index>-1) {
//            cout << "matching: "<<vertexid<<": " << pathid1_index << "," << pathid2_index << endl;
//        }
        if (pathid1_index > -1 && pathid2_index > -1) {
            matchedvertexids.push_back(vertexid);
        }
    }

    return matchedvertexids;

}


vector<vector<int>> AssociationRules::ParallelFindTarget(OriginalGraph & graph, vector<pair<int,int>>& vertexidset, int edgelabel, int length){

    vector<int> temptargets;
    vector<vector<int>> targetslist;
    vector<bool> insertedtemptaregets;
    insertedtemptaregets.resize(graph.number_of_vertices);
    targetslist.clear();


    for (auto vertexid: vertexidset) {

        for (int i = 0; i < graph.number_of_vertices; i++)insertedtemptaregets[i] = false;
        temptargets.clear();

        for (auto target: matchedvertices[vertexid.first].targets[length - 1][vertexid.second]) {// the end vertices of path A from the vertex

            for (auto edgeid: graph.vertex_connecting_edge_list[target]) {// connecting edges to the end vertices

                Edge tempedge = graph.edge_list[edgeid];
                if (tempedge.edgelabel == edgelabel && !insertedtemptaregets[tempedge.dst]) {// the edges has the edge legbel
                    temptargets.push_back(tempedge.dst);
                    insertedtemptaregets[tempedge.dst] = true;
                }
            }
        }

        targetslist.push_back(temptargets);
    }

    return targetslist;

}

bool AssociationRules::CandidatePathMatching(Path& path, vector<int>& pathidentifiers, vector<int>& mergetuples, int position, int start, int end, int length){


    Path tmppath = path;

    for (int for_att = 0; for_att < path.vertextuples[mergetuples[position]].size(); for_att++) {

        Path modifiedpath = tmppath;

        modifiedpath.vertextuples[mergetuples[position]].erase(modifiedpath.vertextuples[mergetuples[position]].begin() + for_att);

        //cout<<"CandidatePathMatching ;"<<position<<" ";
        //modifiedpath.Show();

        if(position==mergetuples.size()-1){
            bool containFlag = false;

            for (int for_path = start; for_path < end; for_path++) {

                if (paths[length][for_path] == modifiedpath) {
                    //paths[length][for_path].Show();
                    pathidentifiers.push_back(paths[length][for_path].pathidentifier);
                    //cout<<"CandidatePathMatching ;"<<paths[length][for_path].pathidentifier<<endl;
                    containFlag = true;
                    break;
                }
            }
            if (!containFlag) {
                pathidentifiers.clear();
                return false;
            }
        }
        else {
            int nextposition = position+1;
            if(!CandidatePathMatching(modifiedpath, pathidentifiers, mergetuples, nextposition, start, end, length)) return false;
        }
    }

    return true;
}


vector<pair<int,int>> AssociationRules::ParallelPathContainment(vector<int>& vertexidset, int length, int pathidentifier){

    int pathid;
    vector<pair<int,int>> matchedvertexsets;
    matchedvertexsets.clear();

    for (int vertexid: vertexidset) {
        pathid = Contains_ReturnIndex(matchedvertices[vertexid].pathidentifiers[length], pathidentifier);//vertex matching the path A
        if (pathid > -1)matchedvertexsets.push_back(make_pair(vertexid,pathid));
    }

    return matchedvertexsets;
}


void AssociationRules::VerticesPartitioning(OriginalGraph & graph, int mode, bool initialFlag){

    if(initialFlag) {
        parallelvertexidsets.clear();
        parallelvertexidsets.resize(corenum);
        for (int i = 0; i < graph.number_of_vertices; i++) {
            parallelvertexidsets[i % corenum].push_back(i);
        }
    }
    else if(mode==RANDOM) {
        parallelvertexidsets.clear();
        parallelvertexidsets.resize(corenum);
        int vertexcount=0;
        if(StratifiedSmaplingFlag) {
            for (int for_stratum = 0; for_stratum < graph.strata.size(); for_stratum++) {

                bool candidateFlag = false;
                for (auto attribute: frequentattributes) {
                    if (Contains(graph.strataid2attributes[for_stratum], attribute)) {
                        candidateFlag = true;
                        break;
                    }
                }

                if (candidateFlag) {

                    for (auto vertexid: graph.strata[for_stratum]) {
                        if (rand() / (double) RAND_MAX > ApproximateSamplingRate)continue;
                        parallelvertexidsets[vertexcount % corenum].push_back(vertexid);
                        vertexcount++;
                    }

                }

            }
        }
        else {
            for (int i = 0; i < graph.number_of_vertices; i++) {
                if (rand() / (double) RAND_MAX > ApproximateSamplingRate)continue;
                parallelvertexidsets[vertexcount % corenum].push_back(i);
                vertexcount++;
            }
        }
    }
    else if(mode==AVERAGESCORE){
        vector<double> corescores(corenum, 0);
        vector<double> corevertexnum(corenum, 0);
        vector<pair<double,int>> vertexscores;

        parallelvertexidsets.clear();
        parallelvertexidsets.resize(corenum);

        if(StratifiedSmaplingFlag) {
            for (int for_stratum = 0; for_stratum < graph.strata.size(); for_stratum++) {

                bool candidateFlag = false;
                for (auto attribute: frequentattributes) {
                    if (Contains(graph.strataid2attributes[for_stratum], attribute)) {
                        candidateFlag = true;
                        break;
                    }
                }

                if (candidateFlag) {

                    for (auto vertexid: graph.strata[for_stratum]) {
                        if (rand() / (double) RAND_MAX > ApproximateSamplingRate) {
                            vertexscores.push_back(make_pair(0, vertexid));
                            continue;
                        }

                        double score;
                        score = matchedvertices[vertexid].pathidentifiers[0].size();
                        double edge_score = 0;
                        for (auto &edgelabel: candidateedgelabels[0]) {
                            edge_score += graph.vertex_connectingcount_edgelabels[vertexid][edgelabel];
                        }
                        //cout<<score<<","<<edge_score<<endl;
                        vertexscores.push_back(make_pair(score * edge_score, (int) vertexid));
                    }

                }

            }
        }
        else {
            for (int i = 0; i < graph.number_of_vertices; i++) {
                if (rand() / (double) RAND_MAX > ApproximateSamplingRate) {
                    vertexscores.push_back(make_pair(0, i));
                    continue;
                }

                double score;
                score = matchedvertices[i].pathidentifiers[0].size();
                double edge_score = 0;
                for (auto &edgelabel: candidateedgelabels[0]) {
                    edge_score += graph.vertex_connectingcount_edgelabels[i][edgelabel];
                }
                //cout<<score<<","<<edge_score<<endl;
                vertexscores.push_back(make_pair(score * edge_score, i));
            }
        }

        sort(vertexscores.begin(),vertexscores.end(), std::greater<pair<double,int>>());


        for (auto& vertexscore : vertexscores) {

            if(vertexscore.first==0)break;
            double minscore = graph.number_of_vertices;
            int minindex = 0;

            for(int for_core=0;for_core<corenum;for_core++) {
                if(corescores[for_core]==0){
                    minindex=for_core;
                    break;
                }
                else if(corescores[for_core]<minscore){
                    minindex=for_core;
                    minscore=corescores[for_core];
                }
            }

            //cout<<minscore<<","<<minindex<<endl;
            parallelvertexidsets[minindex].push_back(vertexscore.second);
            corescores[minindex]+=vertexscore.first;
            corevertexnum[minindex]++;
        }
    }
    else if(mode==AVERAGESCOREVERTEXNUMBALANCE){
        vector<double> corescores(corenum, 0);
        vector<double> corevertexnum(corenum, 0);
        vector<pair<double,int>> vertexscores;

        parallelvertexidsets.clear();
        parallelvertexidsets.resize(corenum);

        if(StratifiedSmaplingFlag) {
            for (int for_stratum = 0; for_stratum < graph.strata.size(); for_stratum++) {

                bool candidateFlag = false;
                for (auto attribute: frequentattributes) {
                    if (Contains(graph.strataid2attributes[for_stratum], attribute)) {
                        candidateFlag = true;
                        break;
                    }
                }

                if (candidateFlag) {

                    for (auto vertexid: graph.strata[for_stratum]) {
                        if (rand() / (double) RAND_MAX > ApproximateSamplingRate) {
                            vertexscores.push_back(make_pair(0, vertexid));
                            continue;
                        }

                        double score;
                        score = matchedvertices[vertexid].pathidentifiers[0].size();
                        double edge_score = 0;
                        for (auto &edgelabel: candidateedgelabels[0]) {
                            edge_score += graph.vertex_connectingcount_edgelabels[vertexid][edgelabel];
                        }
                        //cout<<score<<","<<edge_score<<endl;
                        vertexscores.push_back(make_pair(score * edge_score, (int) vertexid));
                    }

                }

            }
        }
        else {
            for (int i = 0; i < graph.number_of_vertices; i++) {

                if (rand() / (double) RAND_MAX > ApproximateSamplingRate) {
                    vertexscores.push_back(make_pair(0, i));
                    continue;
                }
                double score;
                score = matchedvertices[i].pathidentifiers[0].size();
                double edge_score = 0;
                for (auto &edgelabel: candidateedgelabels[0]) {
                    edge_score += graph.vertex_connectingcount_edgelabels[i][edgelabel];
                }
                vertexscores.push_back(make_pair(score * edge_score, i));
            }
        }

        sort(vertexscores.begin(),vertexscores.end(), std::greater<pair<double,int>>());

        for (auto& vertexscore : vertexscores) {

            if(vertexscore.first==0)break;
            double minscore = graph.number_of_vertices;
            int minscoreindex = 0;
            int minvertexnum=graph.number_of_vertices;
            int maxvertexnum=0;
            int minvertexindex=0;

            for(int for_core=0;for_core<corenum;for_core++) {
                if(corevertexnum[for_core]==0){
                    minscoreindex=for_core;
                    minvertexnum=0;
                    break;
                }
                else if(corescores[for_core]<minscore){
                    minscoreindex=for_core;
                    minscore=corescores[for_core];
                }

                if(minvertexnum>corevertexnum[for_core]){
                    minvertexnum=corevertexnum[for_core];
                    minvertexindex=for_core;
                }
                else if(maxvertexnum<corevertexnum[for_core]){
                    maxvertexnum=corevertexnum[for_core];
                }
            }
            if(maxvertexnum-minvertexnum > 2){
                parallelvertexidsets[minvertexindex].push_back(vertexscore.second);
                corescores[minvertexindex]+=vertexscore.first;
                corevertexnum[minvertexindex]++;

            }
            else {
                parallelvertexidsets[minscoreindex].push_back(vertexscore.second);
                corescores[minscoreindex]+=vertexscore.first;
                corevertexnum[minscoreindex]++;

            }
        }
    }


    int vertexcountcheck=0;
    for(int for_core=0;for_core<corenum;for_core++) {
        cout<<"core"<<for_core<<": "<<parallelvertexidsets[for_core].size()<<endl;
    }
//    if(vertexcountcheck!=graph.number_of_vertices){
//        cout<<"partition error!!"<<endl;
//    }
}


void AssociationRules::MetricsComputation(OriginalGraph& graph){


    Pattern pattern;
    int pathid1;
    int pathid2;
     pair<int, int> pathindex1;
     pair<int, int> pathindex2;
     Path path1;
     Path path2;


    int patternsize=patterns.size();
    #pragma omp parallel for private(pattern, pathid1, pathid2, pathindex1, pathindex2, path1, path2)
    for(int for_pattern=0;for_pattern<patternsize;for_pattern++){
    //for(auto& pattern: patterns){
        pattern =patterns[for_pattern];
        //if(pattern.duplicateFlag)continue;

        pathid1=pattern.pathidentifiers[0];
        pathid2=pattern.pathidentifiers[1];

        pathindex1 = pathidentifier2pathsindex[pathid1];
        pathindex2 = pathidentifier2pathsindex[pathid2];

        path1 = paths[pathindex1.first][pathindex1.second];
        path2 = paths[pathindex2.first][pathindex2.second];

        //cout<<" asup="<<pattern.count<<", count1="<<path1.count<<", count2="<<path2.count<<endl;

        patterns[for_pattern].confidence1=(double)pattern.count/(double)path1.count;
        patterns[for_pattern].confidence2=(double)pattern.count/(double)path2.count;

        patterns[for_pattern].lift =((double)pattern.count*(double)graph.number_of_vertices)/((double)path1.count*(double)path2.count);

        //cout<<" asup="<<pattern.count<<", conf1="<< pattern.confidence1<<", conf2="<< pattern.confidence2<<endl;

    }

}

void AssociationRules::DuplicateRemove(){


//    for(int i=0;i<patterns.size();i++) {
//        if(patterns[i].pathidentifiers[0]==patterns[i].pathidentifiers[1]){
//            patterns.erase(patterns.begin()+i);
//            i--;
//        }
//    }


    Path path1;
    Path path2;
    int j;

    #pragma omp parallel for private(j)
    for(int i=0;i<patterns.size();i++) {
        if(patterns[i].duplicateFlag)continue;
//
//        path1 = paths[pathidentifier2pathsindex[patterns[i].pathidentifiers[0]].first][pathidentifier2pathsindex[patterns[i].pathidentifiers[0]].second];
//        path2 = paths[pathidentifier2pathsindex[patterns[i].pathidentifiers[1]].first][pathidentifier2pathsindex[patterns[i].pathidentifiers[1]].second];
//
//        if(path2>path1){
//            swap(patterns[i].pathidentifiers[0],patterns[i].pathidentifiers[1]);
//        }

        for (j = i + 1; j < patterns.size(); j++) {

            if(patterns[i]==patterns[j]){
                patterns[j].duplicateFlag=true;
            }
        }
    }

}



void AssociationRules::ShowAllPaths(string filename){


    //string filename="./result/patterns_"+dataname+"_"+to_string(minsupport)+"_"+to_string(maximumpathlength);
    std::ofstream fout(filename, ios::app);

    fout<<"All paths output"<<endl;
    for(int for_length=0;for_length<paths.size();for_length++){

        for(auto path: paths[for_length]){

            fout<<path.pathidentifier<<": ";
            for(int for_vertex=0;for_vertex<for_length+1;for_vertex++){

                fout<<"( ";
                for(int for_tuple=0;for_tuple<path.vertextuples[for_vertex].size();for_tuple++){

                    fout<<path.vertextuples[for_vertex][for_tuple]<<" ";

                }
                fout<<") ";

                if(for_vertex!=for_length){
                    if(path.kleenestar)fout<<path.edgelabels[for_vertex]<<"* ";
                    else fout<<path.edgelabels[for_vertex]<<" ";
                }
            }

            fout<<"count="<<path.count<<endl;
        }
    }

}


void AssociationRules::ShowAllPatterns(string filename){

    //string filename="./result/patterns_"+dataname+"_"+to_string(minsupport)+"_"+to_string(maximumpathlength);
    std::ofstream fout(filename, ios::app);

    fout<<"All pattern output"<<endl;
    for(auto pattern: patterns){
        if(pattern.duplicateFlag)continue;
        fout<<pattern.patternidentifier<<": ";

        Path path = paths[pathidentifier2pathsindex[pattern.pathidentifiers[0]].first][pathidentifier2pathsindex[pattern.pathidentifiers[0]].second];
        fout<<path.pathidentifier<<". ";


        for(int for_vertex=0;for_vertex<path.vertextuples.size();for_vertex++){

            fout<<"( ";
            for(int for_tuple=0;for_tuple<path.vertextuples[for_vertex].size();for_tuple++){

                fout<<path.vertextuples[for_vertex][for_tuple]<<" ";

            }
            fout<<") ";

            if(for_vertex!=path.vertextuples.size()-1){
                if(path.kleenestar)fout<<path.edgelabels[for_vertex]<<"* ";
                else fout<<path.edgelabels[for_vertex]<<" ";
            }
        }
        path = paths[pathidentifier2pathsindex[pattern.pathidentifiers[1]].first][pathidentifier2pathsindex[pattern.pathidentifiers[1]].second];
        fout<<" , "<<path.pathidentifier<<". ";


        for(int for_vertex=0;for_vertex<path.vertextuples.size();for_vertex++){

            fout<<"( ";
            for(int for_tuple=0;for_tuple<path.vertextuples[for_vertex].size();for_tuple++){

                fout<<path.vertextuples[for_vertex][for_tuple]<<" ";

            }
            fout<<") ";

            if(for_vertex!=path.vertextuples.size()-1){
                if(path.kleenestar)fout<<path.edgelabels[for_vertex]<<"* ";
                else fout<<path.edgelabels[for_vertex]<<" ";
            }
        }


        fout<<" asup="<<pattern.count<<", conf1="<<pattern.confidence1<<", conf2="<<pattern.confidence2<<", lift="<<pattern.lift<<endl;
    }


}



void AssociationRules::SimpleShowAllPatterns(string filename){

    //string filename="./result/patterns_"+dataname+"_"+to_string(minsupport)+"_"+to_string(maximumpathlength);
    std::ofstream fout(filename, ios::app);

    fout<<"All pattern output"<<endl;
    for(auto pattern: patterns){
        if(pattern.duplicateFlag)continue;

        Path path = paths[pathidentifier2pathsindex[pattern.pathidentifiers[0]].first][pathidentifier2pathsindex[pattern.pathidentifiers[0]].second];

        for(int for_vertex=0;for_vertex<path.vertextuples.size();for_vertex++){

            fout<<"( ";
            for(int for_tuple=0;for_tuple<path.vertextuples[for_vertex].size();for_tuple++){

                fout<<path.vertextuples[for_vertex][for_tuple]<<" ";

            }
            fout<<") ";

            if(for_vertex!=path.vertextuples.size()-1){
                if(path.kleenestar)fout<<path.edgelabels[for_vertex]<<"* ";
                else fout<<path.edgelabels[for_vertex]<<" ";
            }
        }
        path = paths[pathidentifier2pathsindex[pattern.pathidentifiers[1]].first][pathidentifier2pathsindex[pattern.pathidentifiers[1]].second];
        fout<<"\t";


        for(int for_vertex=0;for_vertex<path.vertextuples.size();for_vertex++){

            fout<<"( ";
            for(int for_tuple=0;for_tuple<path.vertextuples[for_vertex].size();for_tuple++){

                fout<<path.vertextuples[for_vertex][for_tuple]<<" ";

            }
            fout<<") ";

            if(for_vertex!=path.vertextuples.size()-1){
                if(path.kleenestar)fout<<path.edgelabels[for_vertex]<<"* ";
                else fout<<path.edgelabels[for_vertex]<<" ";
            }
        }

        fout<<endl;
    }

}

// combine two frequent paths to obtain longer paths: incorrect in our definition
//void AssociationRules::FindFrequentLengthLPath_test(OriginalGraph& graph, int length) {
//
//
//    int countmax=0;
//
//    countmax=0;
//    //cout<<"START finding length "<< length << " frequent paths"<<endl;
//
//
//    for(auto pathL: paths[length-1]) {
//
//        for (auto path_1: paths[1]) {
//
//            if (pathL.vertextuples.back() != path_1.vertextuples[0])continue;
//
//            //cout<<"("<<pathL.vertextuples[0][0]<<") "<< pathL.edgelabels[0] <<" ( "<<pathL.vertextuples[1][0]<<" ) "<< path_1.edgelabels[0]<<" ( " <<path_1.vertextuples[1][0]<<")"<<endl;
//
//            int count = 0;
//            vector<int> matchedvertexids;
//            matchedvertexids.clear();
//
//            vector<pair<int, int>> matchedvertex_pathidindexs;
//            matchedvertexids.clear();
//            for (auto matchedvertex: matchedvertices) {//for all vertices
//                int pathid = Contains_ReturnIndex(matchedvertex.pathidentifiers[length - 1], pathL.pathidentifier);//vertex matching the path A
//                if (pathid > -1)matchedvertex_pathidindexs.push_back(make_pair(matchedvertex.vertexid, pathid));
//            }
//
//            int numberOfmatchedvertex = matchedvertex_pathidindexs.size();
//            vector<bool> alltargets;
//            alltargets.resize(graph.number_of_vertices);
//            vector<bool> insertedtemptaregets;
//            insertedtemptaregets.resize(graph.number_of_vertices);
//
//            vector<int> temptargets;
//            vector<vector<int>>targetslist;
//
//            targetslist.clear();
//
//            for(int i=0;i<graph.number_of_vertices;i++){
//                alltargets[i]=false;
//                insertedtemptaregets[i]=false;
//            }
//
//            vector<vector<int>> nexttargets;
//            nexttargets.resize(graph.number_of_vertices);
//            vector<bool> matchedvertexresult;
//            matchedvertexresult.resize(graph.number_of_vertices);
//
//            vector<bool> insertedmatchedvertexs;
//            insertedmatchedvertexs.resize(graph.number_of_vertices);
//            count=0;
//            for(int i=0;i<graph.number_of_vertices;i++){
//                nexttargets[i].clear();
//                insertedmatchedvertexs[i]=false;
//                matchedvertexresult[i]=false;
//            }
//
//            for (auto matchedvertex_pathidindex: matchedvertex_pathidindexs){
//                temptargets.clear();
//
//                for (auto target: matchedvertices[matchedvertex_pathidindex.first].targets[length - 1][matchedvertex_pathidindex.second]) {// the end vertices of path A from the vertex
//
//                    int pathid = Contains_ReturnIndex(matchedvertices[target].pathidentifiers[length - 1], path_1.pathidentifier);//vertex matching the path A
//                    if (pathid > -1) {
//                        //if (matchedattributesets_vertex[target][path_1.pathidentifier]) {
//
//                        if (!insertedmatchedvertexs[matchedvertex_pathidindex.first]) {
//                            matchedvertexids.push_back(matchedvertex_pathidindex.first);
//                            count++;
//                            insertedmatchedvertexs[matchedvertex_pathidindex.first] = true;
//                        }
//                        for (int i = 0; i < matchedvertices[target].targets[1][pathid].size(); i++) {
//                            nexttargets[matchedvertex_pathidindex.first].push_back(matchedvertices[target].targets[1][pathid][i]);
//                        }
//                    }
//                }
//            }
//
//            for(auto& targets: nexttargets) {
//                std::sort(targets.begin(), targets.end());
//                targets.erase(std::unique(targets.begin(), targets.end()), targets.end());
//            }
////
////            for(int i=0;i<graph.number_of_vertices;i++){
////                nexttargets[i].clear();
////                insertedmatchedvertexs[i]=false;
////                matchedvertexresult[i]=false;
////            }
////
////
////            for (int for_vertex=0; for_vertex < numberOfmatchedvertex; for_vertex++) {
////
////                for (auto target: targetslist[for_vertex]) {
////
////                    if (matchedattributesets_vertex[target][path_1.pathidentifier]) {
////
////                        if (!insertedmatchedvertexs[matchedvertex_pathidindexs[for_vertex].first]) {
////                            matchedvertexids.push_back(matchedvertex_pathidindexs[for_vertex].first);
////                            count++;
////                            insertedmatchedvertexs[matchedvertex_pathidindexs[for_vertex].first] = true;
////                        }
////                        nexttargets[matchedvertex_pathidindexs[for_vertex].first].push_back(target);
////                    }
////
////                }
////            }
//
//            //cout<<"count="<<count<<endl;
//            if(count>countmax)countmax=count;
//
//            if (count > minsupport) {
//
//                Path temppath;
//                temppath.pathidentifier = globalpathidcount;
//                globalpathidcount++;
//                temppath.count = count;
//                temppath.vertextuples = pathL.vertextuples;
//                temppath.vertextuples.push_back(path_1.vertextuples[1]);
//                temppath.edgelabels = pathL.edgelabels;
//                temppath.edgelabels.push_back(path_1.edgelabels[0]);
//                temppath.kleenestar = false;
//
//                paths[length].push_back(temppath);
//
//                for (auto matchedvertexid: matchedvertexids) {
//
//                    matchedvertices[matchedvertexid].pathidentifiers[length].push_back(temppath.pathidentifier);
//                    matchedvertices[matchedvertexid].targets[length].push_back(nexttargets[matchedvertexid]);
//                }
//            }
//        }
//    }
//
//}