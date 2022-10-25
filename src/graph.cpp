
#include "graph.h"

void OriginalGraph::InputFile(std::string graph_directory_name)
{
    std::string vertexfile_name= graph_directory_name+"/vertex";
    std::string edgefile_name= graph_directory_name+"/edge";

    const char *vertexfile=vertexfile_name.c_str();
    const char *edgefile=edgefile_name.c_str();
    std::ifstream inputvertex(vertexfile);

    if(!inputvertex){
        std::cout<<"error: cannot open vertex file:"+vertexfile_name <<std::endl;
        exit(1);
    }

    std::ifstream inputedge(edgefile);

    if(!inputedge){
        std::cout<<"error: cannot open edge files:"+edgefile_name<<std::endl;
        exit(1);
    }

    //inputvertex>>number_of_vertices>>cardinality_of_vertextuples;
    //cardinality_of_vertextuples.resize(number_of_tuples);

    number_of_vertices=0;
    number_of_edges=0;
    cardinality_of_edgelabels=0;
    cardinality_of_vertextuples=0;
    average_attribute=0;

    string tmp;
    string line;
    char del = ' ';
    bool firstline=true;
    vertexname2vertexid.set_empty_key("-1");
    vertexid2vertexname.set_empty_key(-1);
	attribute2attributeid.set_empty_key("-1");
	attributeid2attribute.set_empty_key(-1);

    cout<<"vertex input START"<<endl;
    while(getline(inputvertex, line)){
        if(line.empty())continue;

        int tmp_vertexid;
        vector<int> tuple;
        bool firstelement=true;
        for (const auto element : split(line, del)) {

            //cout<<line<<","<<element<<endl;
            if(firstelement){
                vertexname2vertexid.insert({element,number_of_vertices});
                vertexid2vertexname.insert({number_of_vertices,element});
                tmp_vertexid=vertexname2vertexid[element];

                number_of_vertices++;
                firstelement=false;
            }
            else {
                if(attribute2attributeid.find(element)==attribute2attributeid.end()){
                    attribute2attributeid.insert({element,cardinality_of_vertextuples});
                    attributeid2attribute.insert({cardinality_of_vertextuples,element});
                    cardinality_of_vertextuples++;
                }

                tuple.push_back(attribute2attributeid[element]);
                average_attribute++;
            }
        }
        sort(tuple.begin(),tuple.end());
        Vertex vertex(tmp_vertexid, tuple);
        vertex_list.push_back(vertex);

        bool sameFlag=false;
        for(auto att2strata: attributes2strataid){
            if(tuple.size()!=att2strata.first.size())continue;
            sameFlag=true;
            for(int for_att=0;for_att<att2strata.first.size();for_att++){
                if(tuple[for_att]!=att2strata.first[for_att]){
                    sameFlag=false;
                    break;
                }
            }
            if(sameFlag){
                strata[att2strata.second].push_back(tmp_vertexid);
                break;
            }
        }
        if(!sameFlag){
            attributes2strataid.push_back(make_pair(tuple,strata.size()));

            strata.push_back({tmp_vertexid});
            strataid2attributes.push_back(tuple);
        }



    }


    numbers_of_attributes.resize(cardinality_of_vertextuples);
    for(int i=0;i<cardinality_of_vertextuples;i++)numbers_of_attributes[i]=0;
    for(auto& vertex: vertex_list){
        for(auto& attribute: vertex.vertextuple){
            numbers_of_attributes[attribute]++;
        }
    }



    string v1,v2,label;
    //inputedge>>number_of_edges>>cardinality_of_edgelabels;
    vertex_connecting_edge_list.resize(number_of_vertices);
    vertex_inverseconnecting_edge_list.resize(number_of_vertices);
    edgelabel2edgelabelid.set_empty_key("-1");
	edgelabelid2edgelabel.set_empty_key(-1);

    for(int i=0;i<number_of_vertices;i++){
        vertex_connecting_edge_list[i].clear();
        vertex_inverseconnecting_edge_list[i].clear();
    }

    cout<<"vertex input DONE"<<endl;

    cout<<"edge input START"<<endl;
    while(!inputedge.eof()){

        inputedge >> v1 >> label >> v2;
        if(inputedge.fail())break;

        if(edgelabel2edgelabelid.find(label)==edgelabel2edgelabelid.end()){
            edgelabel2edgelabelid.insert({label,cardinality_of_edgelabels});
            edgelabelid2edgelabel.insert({cardinality_of_edgelabels,label});
            cardinality_of_edgelabels++;
        }


        Edge edge(edge_list.size(), vertexname2vertexid[v1], vertexname2vertexid[v2], edgelabel2edgelabelid[label]);


        edge_list.push_back(edge);
        vertex_connecting_edge_list[vertexname2vertexid[v1]].push_back(edge.edgeid);
        vertex_inverseconnecting_edge_list[vertexname2vertexid[v2]].push_back(edge.edgeid);
        number_of_edges++;
    }

    cout<<"edge input DONE"<<endl;

    cout<<"statics compute START"<<endl;
    cout<<"cardinality of edge lages = "<<cardinality_of_edgelabels<<", cardinality of vertex attributes = "<<cardinality_of_vertextuples<<endl;
    numbers_of_edgelabels.resize(cardinality_of_edgelabels);
    for(int i=0;i<cardinality_of_edgelabels;i++)numbers_of_edgelabels[i]=0;

    vertex_connectingcount_edgelabels.resize(number_of_vertices);
    vertex_inverseconnectingcount_edgelabels.resize(number_of_vertices);

    for(int i=0;i<number_of_vertices;i++){
        vertex_connectingcount_edgelabels[i].resize(cardinality_of_edgelabels);
        vertex_inverseconnectingcount_edgelabels[i].resize(cardinality_of_edgelabels);
        for(int j=0;j<cardinality_of_edgelabels;j++){
            vertex_connectingcount_edgelabels[i][j]=0;
            vertex_inverseconnectingcount_edgelabels[i][j]=0;
        }
    }

    //cout<<"statics2 compute DONE"<<endl;
    attribute_connecting_vertex.resize(cardinality_of_vertextuples);
    attribute_inverseconnectingcount_edgelabels.resize(cardinality_of_vertextuples);
    attribute_max_inverseconnectingcount_edgelabels.resize(cardinality_of_vertextuples);


    for(int i=0;i < cardinality_of_vertextuples;i++){
        attribute_connecting_vertex[i].clear();
        attribute_inverseconnectingcount_edgelabels[i].resize(cardinality_of_edgelabels);
        attribute_max_inverseconnectingcount_edgelabels[i]=0;
        for(int j=0;j < cardinality_of_edgelabels;j++){
            attribute_inverseconnectingcount_edgelabels[i][j]=0;
        }
    }
    //cout<<"statics3 compute DONE"<<endl;

    for(auto& edge: edge_list){
        numbers_of_edgelabels[edge.edgelabel]++;
        vertex_connectingcount_edgelabels[edge.src][edge.edgelabel]++;
        vertex_inverseconnectingcount_edgelabels[edge.dst][edge.edgelabel]++;

        for(auto& attribute: vertex_list[edge.dst].vertextuple) {

            attribute_inverseconnectingcount_edgelabels[attribute][edge.edgelabel]++;
        }
    }

    //cout<<"statics4 compute DONE"<<endl;
    for(int i=0;i < cardinality_of_vertextuples;i++){
        for(int j=0;j < cardinality_of_edgelabels;j++){
            if(attribute_max_inverseconnectingcount_edgelabels[i] < attribute_inverseconnectingcount_edgelabels[i][j]){
                attribute_max_inverseconnectingcount_edgelabels[i] = attribute_inverseconnectingcount_edgelabels[i][j];
            }
        }
    }

    connecting_vertex.clear();
    for(auto& vertex : vertex_list){
        for(auto& attribute : vertex.vertextuple){
            attribute_connecting_vertex[attribute].push_back(vertex.vertexid);
        }
        connecting_vertex.push_back(vertex.vertexid);
    }

    //cout<<"statics5 compute DONE"<<endl;

    maximumIndegree_edgelabel=0;
    averageIndegree=0;
    double count_for=0;
    for(auto& vertex_inverseconnectingcount: vertex_inverseconnectingcount_edgelabels){

        for(auto& count: vertex_inverseconnectingcount){

            if(maximumIndegree_edgelabel<count)maximumIndegree_edgelabel=count;
            averageIndegree+=count;
            count_for++;
        }
    }
    averageIndegree= averageIndegree / number_of_vertices;
    average_attribute=average_attribute/(double)number_of_vertices;
    //cout<<maximumIndegree_edgelabel<<endl;




}