#include <iostream>

#include "../DynamicCover.cpp"
#include "../DynamicGraph.cpp"

#include <NetworKit/io/EdgeListWriter.h>
#include <NetworKit/io/CoverWriter.h>

using namespace std;

DynamicGraph dg;
DynamicCover dc;

int T = 1000;

void readDynamicGraph(){
	dg.readDynamicGraph("ckbDynamicGraphByNode");
}

void readDynamicCover(){
	dc.readCover("ckbDynamicCoverByCommunity");
}


int main(){
	readDynamicGraph();
	readDynamicCover();
	int t = 0;
	NetworKit::Graph g = dg.snapshot(t);
	NetworKit::Cover c0 = dc.snapshot(t, g.numberOfNodes());
	while (t <= T){
		g = dg.snapshot(t);
		NetworKit::Cover c = dc.snapshot(t, g.numberOfNodes());
		auto x = c.subsetSizes();
		sort(x.begin(), x.end());
		int min = x[0], max = x[x.size()-1];
		int med = x[x.size()/2];
		int avg = 0;
		for (int k = 0; k < x.size(); k++) avg += x[k];
		avg /= x.size();

		cout << t << " "  << min << " " << max << " " << med << " " << avg << " " << x.size() << endl;
		t += 1;
	}
}
