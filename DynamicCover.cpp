#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <string.h>
#include <cstring>

#include <NetworKit/io/EdgeListReader.h>
#include <NetworKit/io/EdgeListWriter.h>
#include <NetworKit/io/EdgeListPartitionReader.h>
#include <NetworKit/structures/Cover.h>

#include <NetworKit/io/CoverReader.h>

using namespace std;

class DynamicCover{
private:
	class CommunityAssign{
	public:
		NetworKit::index communityId;
		int start;
		int end;
	};

public:
	map<NetworKit::node, vector<CommunityAssign>> dynCover;

	void readCover(string coverFileName){
		ifstream coverFile(coverFileName);
		NetworKit::index commId;
		string line;
		while (getline(coverFile, line)){
			if (line.find(':') != string::npos){
				char *cstr = &line[0u];
				char *t = strtok(cstr,":");
				commId = atoi(t);

				cstr = strtok(NULL, ":");
				t = strtok(cstr, "(");
				t = strtok(t, ",");
				NetworKit::node nodeId = atoi(t);
				CommunityAssign ca;
				ca.communityId = commId;
				t = strtok(NULL, ",");
				ca.start = atoi(t);
				t = strtok(NULL,")");
				ca.end = atoi(t);
				dynCover[nodeId].push_back(ca);
			}
			else if (line.find(',') != string::npos){
				char *cstr = &line[0u];
				char *t = strtok(cstr, "(");
				t = strtok(t, ",");
				CommunityAssign ca;
				NetworKit::node nodeId = atoi(t);
				ca.communityId = commId;
				t = strtok(NULL, ",");
				ca.start = atoi(t);
				t = strtok(NULL,")");
				ca.end = atoi(t);
				dynCover[nodeId].push_back(ca);

			}
		}
	}

	NetworKit::Cover snapshot(int timestamp, NetworKit::index numElements){
		NetworKit::Cover cover(numElements);
		for (auto it = dynCover.begin(); it != dynCover.end(); ++it){
			NetworKit::node nodeId = it->first;
			auto comms = it->second;
			for (auto it1 = comms.begin(); it1 != comms.end(); ++it1){
				CommunityAssign ca = *it1;
				if (ca.start <= timestamp && ca.end >= timestamp){
					NetworKit::index commId = ca.communityId;
					if (cover.upperBound() <= commId)
						cover.setUpperBound(commId + 1);
					cover.addToSubset(commId, nodeId);
					//cout << " Adding " << nodeId << " to " << commId << endl;
				}
			}
		}
		return cover;
	}

	double computeF1(NetworKit::Graph graph, NetworKit::Cover s1, NetworKit::Cover s2, bool debug){
		/*Construct GroundTruth sets*/
		vector<vector<NetworKit::node> > s1Sets (s1.upperBound());
		vector<NetworKit::count> s1Sizes(s1.upperBound(), 0), s2Sizes(s2.upperBound(), 0);

		graph.forNodes([&](NetworKit::node u) {
			for (NetworKit::index s : s1[u]) {
				s1Sets[s].push_back(u);
				++s1Sizes[s];
			}
			for (NetworKit::index s : s2[u]) {
				++s2Sizes[s];
			}
		});
		map<NetworKit::index, NetworKit::count> overlap;
		map<NetworKit::index, double> f1Agreement;

		for (NetworKit::index gt = 0; gt < s1Sets.size(); ++gt){
			if (s1Sets[gt].size() > 0) {
				overlap.clear();
				auto &s1Nodes = s1Sets[gt];

				for (NetworKit::node u : s1Nodes)
					for (NetworKit::index s : s2[u])
						++overlap[s];
				double bestF1 = 0;
				int bestSize = 0;
				double bestP = 0, bestR = 0;
				for (auto o : overlap) {
					double precision = o.second * 1.0 / s2Sizes[o.first];
					double recall = o.second * 1.0 / s1Nodes.size();
					double f1 = 2 * (precision * recall) / (precision + recall);
					if (f1 > bestF1){
						bestF1 = f1;
						bestSize = s2Sizes[o.first];
						bestP = precision;
						bestR = recall;
					}
				}
				if (debug){
				//	cout << "[f1] OUT:\t" << s1Nodes.size() << "\t" << bestSize << "\t" << bestP << "\t" << bestR << "\t" << bestF1 << endl;
					if (s1Nodes.size() == 2) cout << "Comm: " << s1Nodes[0] << ", " << s1Nodes[1] << endl;
				}
				//else cout << "[f1] OUT1:\t" << s1Nodes.size() << "\t" << bestSize << "\t" << bestP << "\t" << bestR << "\t" << bestF1 << endl;

				f1Agreement[gt] = bestF1;
			}
		}

		double sumF1 = 0;
		double weightedAverageF1 = 0;
		NetworKit::count sumTruth = 0;

		for (auto a : f1Agreement) {
			sumF1 += a.second;
			if (s1Sizes[a.first] <= 3) continue;
			weightedAverageF1 += a.second * s1Sizes[a.first];
			sumTruth += s1Sizes[a.first];
		}

		weightedAverageF1 = weightedAverageF1 / sumTruth;
		return weightedAverageF1;

	}

	double measureF1AtTimestamp(NetworKit::Graph graph, NetworKit::Cover foundCover, int timestamp){
		NetworKit::Cover trueCover = snapshot(timestamp, graph.numberOfNodes());
		//cout << "f1 :: Number of subsets in trueCover, and foundCover at " << timestamp << " = " << trueCover.numberOfSubsets() << ", " << foundCover.numberOfSubsets() << endl;
		double f1 = computeF1(graph, foundCover, trueCover, true);
		//cout << "I. f1 = " << f1 << endl;
		f1 += computeF1(graph, trueCover, foundCover, false);
		//cout << "II. f1 = " << f1 << endl;

		return 0.5*f1;
	}

};


/*
int main(){
	DynamicCover dc;
	dc.readCover("ckbDynamicCommunitiesO");
	NetworKit::EdgeListReader reader = NetworKit::EdgeListReader('\t',0);
        NetworKit::Graph graph = reader.read("ckbEdgeList.small");

	NetworKit::CoverReader cr = NetworKit::CoverReader();
        NetworKit::Cover trueCover = cr.read("ckbCover.small", graph);

	cout << "F1 (1) = " << dc.measureF1AtTimestamp(graph, trueCover, 0) << endl;
}
*/
