#include "graph.h"
#include "associationrulemining.h"

bool EfficientFlag=true;
bool ParallelFlag=false;
bool ExtentionFlag=true;
bool StratifiedSmaplingFlag=true;
bool EmptysourceFlag=false;
bool NonDominanceOutputFlag=true;
int minimumRulePathLength=1;
double ApproximateCandidateRate=1.0;
double ApproximateSamplingRate=1.0;
int PartitionMode=RANDOM;
double attributetime=0;
double lengthonetime=0;
double lengthLtime=0;
double kleenetime=0;
double combinetime=0;

int main(int argc, char* argv[])
{

   Result result;
    result.clear();
    string resultFile = "./result/test";
    string inputFile;
    int threadnum=1;
    double minimumSupport = 0.1;
    int maximumLength=2;
    double timelimit = 24;
    bool sjoin=true;
    bool stjoin=true;
    bool kleene=true;
    bool patternoutputFlag=true;
    bool attributeoutputFlag=false;
    bool absoluteSupFlag=false;
//    string workloadFile;
//    string graphFile;
//    string indexworkloadFile;
//    string indexFile;
//    int given_k=0;
//    string string_k;
//    bool workloadindex=false;
//    string maintenancegraph;
//    bool maintenancegraphFlag=false;
//    string maintenanceworkload;
//    bool maintenanceworkloadFlag=false;
//    bool maintainedindexoutput=false;

    try {
        for (int i = 1; i < argc; ++i) {
            std::string s = argv[i];
            if (s[0] == '-') {
                s = s.substr(1);

                if (s == "i") {
                    inputFile = argv[++i];
                }
                else if (s == "minsup") {
                    minimumSupport = stod(argv[++i]);
                }
                else if (s == "aminsup") {
                    minimumSupport = stod(argv[++i]);
                    absoluteSupFlag=true;
                }
                else if (s == "o") {
                    resultFile = argv[++i];
                }
                else if (s == "k"){
                    maximumLength = stoi(argv[++i]);
                }
                else if (s == "non"){
                    EfficientFlag=false;
                }
                else if (s == "noextent"){
                    ExtentionFlag=false;
                }
                else if (s == "p") {
                    threadnum = stoi(argv[++i]);
                    ParallelFlag=true;
                }
                else if (s == "pmax") {
                    threadnum = omp_get_max_threads();
                    ParallelFlag=true;
                }
                else if (s == "nooutput") {
                    patternoutputFlag=false;
                }
                else if (s == "attributeoutput") {
                    attributeoutputFlag=true;
                }
                else if (s == "nokleene") {
                    kleene=false;
                }
                else if (s == "limit") {
                    timelimit=stod(argv[++i]);
                }
                else if (s == "cr") {
                    ApproximateCandidateRate=stod(argv[++i]);
                }
                else if (s == "sr") {
                    ApproximateSamplingRate=stod(argv[++i]);
                }
                else if (s == "rs") {
                    StratifiedSmaplingFlag=false;
                }
                else if (s == "es") {
                    EmptysourceFlag=true;
                }
                else if (s == "do") {
                    NonDominanceOutputFlag=false;
                }
                else if (s == "mink"){
                    minimumRulePathLength=stoi(argv[++i]);
                }
                else if (s == "pmode") {
                    string mode=argv[++i];
                    if(mode=="0")PartitionMode=RANDOM;
                    else if(mode=="1")PartitionMode=AVERAGESCORE;
                    else if(mode=="2")PartitionMode=AVERAGESCOREVERTEXNUMBALANCE;
                    else throw std::exception();
                    //PartitionMode=stoi(argv[++i]);
                }
                else {
                    throw std::exception();
                }
            }
            else {
                throw std::exception();
            }
        }
    }
    catch (std::exception& e) {
        cout<<"error input"<<endl;
        return 1;
    }

    vector<string> splitInputFile = split(inputFile, '/');
    string dataname = splitInputFile[splitInputFile.size()-1];

    cout<<"Graph Input Start: "<<inputFile<<endl;
    OriginalGraph originalgraph = OriginalGraph();
    originalgraph.InputFile(inputFile);
    cout<<"Graph Input Done"<<endl;
    if(attributeoutputFlag){
        originalgraph.Outputgraph("./result/"+dataname);
        cout<<"attribute output done"<<endl;
        return 0;
    }


    omp_set_num_threads(threadnum);

    AssociationRules associationrule;
    if(absoluteSupFlag){associationrule.minsupport=minimumSupport;}
    else associationrule.minsupport=minimumSupport*originalgraph.number_of_vertices;

    //associationrule.minsupport=600;
    associationrule.maximumpathlength=maximumLength;
    associationrule.s_join=sjoin;
    associationrule.s_t_join=stjoin;
    associationrule.kleenestar=kleene;
    associationrule.corenum=threadnum;
    associationrule.timelimit=timelimit;

    cout << "Graph statistics: # of vertices " << originalgraph.number_of_vertices << ",  # of edges " << originalgraph.number_of_edges << ", cardinality of edge labels " << originalgraph.cardinality_of_edgelabels << ", cardinality of tuples " << originalgraph.cardinality_of_vertextuples << ", average # of attributes " << originalgraph.average_attribute << ", maximum indegree per edge label " << originalgraph.maximumIndegree_edgelabel << ", average indegree " << originalgraph.averageIndegree << endl;
    cout<<"minsupport = "<<associationrule.minsupport<<" , maxlen = "<<maximumLength<<" , corenum = "<<associationrule.corenum<<endl;

    struct timespec frequentpattern_startTime, frequentpattern_endTime;
    clock_gettime(CLOCK_REALTIME, &frequentpattern_startTime);

    //////Core function/////
    associationrule.FindPattern(originalgraph);
    ////////////////////////

    clock_gettime(CLOCK_REALTIME, &frequentpattern_endTime);
    cout<<"Pattern finding DONE: ";
    double time= (frequentpattern_endTime.tv_sec - frequentpattern_startTime.tv_sec) + (frequentpattern_endTime.tv_nsec - frequentpattern_startTime.tv_nsec) * pow(10, -9);
    cout<<"# of patterns = "<<associationrule.patterns.size()<<", time = "<<time<<endl;

//    struct timespec duplicatecheck_startTime, duplicatecheck_endTime;
//    cout<<"Start Duplicate Check"<<endl;
//    clock_gettime(CLOCK_REALTIME, &duplicatecheck_startTime);
//    associationrule.DuplicateRemove();
//    clock_gettime(CLOCK_REALTIME, &duplicatecheck_endTime);
//    cout<<"DONE: ";
//    cout<<(duplicatecheck_endTime.tv_sec - duplicatecheck_startTime.tv_sec) + (duplicatecheck_endTime.tv_nsec - duplicatecheck_startTime.tv_nsec) * pow(10, -9)<<endl;

    struct timespec confidence_startTime, confidence_endTime;
    cout<<"Start MetricsComputation"<<endl;
    clock_gettime(CLOCK_REALTIME, &confidence_startTime);
    associationrule.MetricsComputation(originalgraph);
    clock_gettime(CLOCK_REALTIME, &confidence_endTime);
    cout<<"DONE: ";
    cout<<(confidence_endTime.tv_sec - confidence_startTime.tv_sec) + (confidence_endTime.tv_nsec - confidence_startTime.tv_nsec) * pow(10, -9)<<endl;



    //cout<<dataname<<endl;
    result.coreNum=associationrule.corenum;
    result.minsup=minimumSupport;
    result.maxlength=maximumLength;
    result.miningTime=time;
    result.attributeTime=attributetime;
    result.lengthoneTime=lengthonetime;
    result.lengthLTime=lengthLtime;
    result.kleeneTime=kleenetime;
    result.combineTime=combinetime;
    result.efficientFlag=EfficientFlag;
    result.extentionFlag=ExtentionFlag;
    result.stratifiedsampleFlag=StratifiedSmaplingFlag;
    result.approximatecandidaterate=ApproximateCandidateRate;
    result.approximatesamplingrate=ApproximateSamplingRate;
    result.inputFile=inputFile;
    result.outputFile=resultFile;
    result.partitionmode=PartitionMode;

    for(int i=0;i < maximumLength+1;i++){
        cout<<associationrule.paths[i].size()<<endl;
        result.pathNum+=associationrule.paths[i].size();
    }

    int patternsize=0;
    for(auto& pattern:associationrule.patterns){
        if(!pattern.duplicateFlag)patternsize++;
    }
    result.patternNum=patternsize;
    cout<<patternsize<<endl;

    result.Output();

    if(patternoutputFlag) {
        string filename="./result/patterns_"+dataname+"_"+to_string(associationrule.minsupport)+"_"+to_string(associationrule.maximumpathlength)+"_"+to_string(result.approximatecandidaterate)+"_"+to_string(result.approximatesamplingrate);

        associationrule.ShowAllPaths(filename);
        associationrule.ShowAllPatterns(filename);

        filename="./result/simplepatterns_"+dataname+"_"+to_string(associationrule.minsupport)+"_"+to_string(associationrule.maximumpathlength)+"_"+to_string(result.approximatecandidaterate)+"_"+to_string(result.approximatesamplingrate);
        associationrule.SimpleShowAllPatterns(filename);
    }


    return 0;
}