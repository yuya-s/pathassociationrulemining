

#ifndef PATHASSOCIATION_GRAPH_H
#define PATHASSOCIATION_GRAPH_H


#include "utility.h"

struct Edge {
public:
    int edgeid;
	int src;
	int dst;
    int edgelabel;

    bool operator==(const Edge& s){return (src == s.src&&dst==s.dst&&edgelabel == s.edgelabel);}
    //bool operator==(const Edge& s){return (edgeid == s.edgeid);}
    //bool operator<(const Edge& s) const {return (src < s.src||(src==s.src&&dst<s.dst)||(src==s.src&&dst==s.dst&&label <= s.label));}

	Edge(int edgeid_,int src_, int dst_,int label_)
	{
        edgeid = edgeid_;
		src = src_;
		dst = dst_;
        edgelabel = label_;
	}

    Edge()
	{
        edgeid = 0;
		src = 0;
		dst = 0;
        edgelabel = 0;
	}
};

struct Vertex {
public:
	int vertexid;
    vector<int> vertextuple;

	Vertex(int vertexid_, vector<int> vertextuple_)
	{
		vertexid = vertexid_;
		vertextuple = vertextuple_;
	}
};


class OriginalGraph {

public:
	int number_of_vertices;
    int number_of_edges;
    //int number_of_tuples;
    int cardinality_of_edgelabels;
    int cardinality_of_vertextuples;
    double average_attribute;

	vector<int> numbers_of_edgelabels;
    vector<int> numbers_of_attributes;
	int maximumIndegree_edgelabel;
    double averageIndegree;


    vector<Vertex> vertex_list;
	vector<Edge> edge_list;
    vector< vector <int> > vertex_connecting_edge_list;
    vector< vector <int> > vertex_inverseconnecting_edge_list;
	vector< vector <int> > vertex_connectingcount_edgelabels;
	vector< vector <int> > vertex_inverseconnectingcount_edgelabels;
    vector< vector <int> > attribute_connecting_vertex;
    vector< int > connecting_vertex; //all vertices ids;
    vector< vector <int> > attribute_inverseconnectingcount_edgelabels;
    vector< int > attribute_max_inverseconnectingcount_edgelabels;

    vector<pair<vector<int>, int>> attributes2strataid;
    vector<vector<int>> strataid2attributes;
    vector<vector<int>> strata; // set of vertex ids

	google::dense_hash_map<std::string, int> vertexname2vertexid;
	google::dense_hash_map<int, std::string> vertexid2vertexname;
	google::dense_hash_map<std::string, int> attribute2attributeid;
	google::dense_hash_map<int, std::string> attributeid2attribute;
	google::dense_hash_map<std::string, int> edgelabel2edgelabelid;
	google::dense_hash_map<int, std::string> edgelabelid2edgelabel;


//    vector< vector <int> > vertex_inverseedge_list;
//	int max_degree;
//    google::dense_hash_map<std::string,int> vertexname2id;
//    google::dense_hash_map<string,int> labelname2id;
//	google::dense_hash_map<int, string> id2vertexname;
//    google::dense_hash_map<int, string> id2labelname;

	OriginalGraph(){
		edge_list.clear();
        vertex_list.clear();
		//max_degree=0;
    }
    void InputFile(std::string);

    void Outputgraph(std::string filename){

        std::ofstream fout(filename, ios::app);
        for(int i=0;i < number_of_vertices;i++){
            fout<<"vid2vname "<<i<<" -> "<<vertexid2vertexname[i]<<endl;
        }
        for(int i=0;i<cardinality_of_edgelabels;i++) {
            fout<<"eid2elabel "<<i<<" -> "<<edgelabelid2edgelabel[i]<<endl;
        }
        for(int i=0;i<cardinality_of_vertextuples;i++) {
            fout<<"aid2attname "<<i<<" -> "<<attributeid2attribute[i]<<endl;
        }
    }
};


#endif //PATHASSOCIATION_GRAPH_H
