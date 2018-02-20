#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <cstring>

#include <NetworKit/graph/Graph.h>

using namespace std;

class DynamicGraph{
private:
	class DynamicEdge{
	public:
		NetworKit::index destination;
		int startTime;
		int endTime;
	};
	class DynamicNode{
	public:
		NetworKit::index nodeId;
		vector<DynamicEdge*> edgeList;
	};
	vector<DynamicNode*> nodeList;
	int maxNodeId;

public:
	void readDynamicGraph(string filename){
		ifstream graphFile(filename);
		NetworKit::index src, dst;
		string line;
		maxNodeId = 0;
		while (getline(graphFile,line)){
			if (line.find(':') != string::npos){
				char *cstr = &line[0u];
				char *t = strtok(cstr,":");
				t = strtok(NULL,":");
				src = atoi(t);
				while (maxNodeId <= src){
					DynamicNode *n = new DynamicNode();
					n->nodeId = maxNodeId;
					nodeList.push_back(n);
					maxNodeId++;
				}
			}
			else{
				char *cstr = &line[0u];
				char *t = strtok(cstr, "(");
				t = strtok(t, ",");
				dst = atoi(t);
				while (maxNodeId <= dst){
					DynamicNode *n = new DynamicNode();
					n->nodeId = maxNodeId;
					nodeList.push_back(n);
					maxNodeId++;
				}
				t = strtok(NULL, ",");
				int st = atoi(t);
				t = strtok(NULL, ",");
				int et = atoi(t);
				t = strtok(NULL,")");
				int commId = atoi(t);
				if (commId == -2 || commId >= 0){
					DynamicEdge* e = new DynamicEdge();
					e->destination = dst;
					e->startTime = st;
					if (et >= 0)
						e->endTime = et;
					else e->endTime = 1001;
					nodeList[src]->edgeList.push_back(e);
				}
			}
		}
		graphFile.close();
	}

	NetworKit::Graph snapshot(int timestamp){
		NetworKit::Graph graph(maxNodeId);
		for (auto n : nodeList){
			for (auto e : n->edgeList){
				if ((e->startTime <= timestamp) && (e->endTime >= timestamp)){
					graph.addEdge(n->nodeId, e->destination);
				}
			}
		}
		return graph;
	}
};

/*
int main(){
	cout << "Hello World!" << endl;
	DynamicGraph dg;
	dg.readDynamicGraph("ckbDynamicGraphByNode");
	NetworKit::Graph g = dg.snapshot(500);

	NetworKit::EdgeListWriter ew(' ', 0);
	ew.write(g, "ckbDynamic500GraphEdgeList");
}
*/
