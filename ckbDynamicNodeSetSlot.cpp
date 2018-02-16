#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <math.h>
#include "PowerlawDegreeSequence.h"
#include "BucketSampling.h"
#include <algorithm>
#include <sstream>

#include <string>
#include <string.h>
/*With node changes*/

using namespace std;

/*************************PARAMETERS******************************/
int T = 50;        //Number of time slots
double lambda = 0.2;//how sharply communities will rise and fall
int N1 = 50000;		//Number of nodes
int N2;             //Number of communities, set in the program
int xmin = 1;		//minimum user-community memberships
int xmax = 50;		//maximum user-community memberships
double beta1 = -2.5;	//community membership exponent
int mmin = 2;		//minimum community size
int mmax = 7500;	    //maximum community size
double beta2 = -2.5;	//community size exponent
double alpha = 2;	//affects intra community edge probability
double gamma_1 = 0.5;	//affects intra community edge probability
double eps = 2;        //inter community edge probability, epsilon = eps/N1
double prob = 4;  //alternative: Intra community edge probability - not used anymore
double probEvent = 0.1; //Probability of an event happening
int timeToMerge = 10;
int timeToSplit = 10;
double minSplitRatio = 0.3; //A community born of split is atleast 0.3 of original community

double probB = 0.25, probD = 0.25, probM = 0.25, probS = 0.25; //individual event probabilities
double initialMem = 0, currMem = 0;
/*******************************************************************/

/*************************DATA STRUCTURES***************************/
struct edge{
	int sourceId;
	int destId;
	int communityId;
	int startTime;
	int endTime;


	edge(int isSourceId, int isDestId, int isCommunityId){
		this->sourceId = isSourceId;
		this->destId = isDestId;
		this->communityId = isCommunityId;
		this->startTime = -1;
		this->endTime = -1;
	}

	void generateStartExp(int startLimit){
		double z = ((double) rand())/((double) RAND_MAX);
		double cumul_prob = 0;
		for (int i = startLimit; i>=0;i--){
			double p_x = lambda * exp( -1* lambda * (startLimit - i));
			cumul_prob += p_x;
			if (z <= cumul_prob){
				this->startTime = i;
				return;
			}
		}
		this->startTime = 0;
	}

	void generateEndExp(int endLimit){
	        if (this->endTime != -1) return;
		double z = ((double) rand())/((double) RAND_MAX);
		double cumul_prob = 0;

		for (int i = endLimit; i<=T;i++){
			double p_x = lambda * exp( -1* lambda * (i - endLimit));
			cumul_prob += p_x;
			if (z <= cumul_prob){
                //cout << "generateEndExp: Setting end time of " << communityId << ": " << sourceId << ", " << destId << ", " << startTime << ", " << endTime << " to " << i << endl << flush;
                if (i < this->startTime) i = this->startTime + 1;
				this->endTime = i;

				return;
			}
		}
		this->endTime = T;
	}
};

struct node{
	int nodeId;
	vector<edge*> adj;
	vector<int> communities;
	int startTime;
	int endTime;
	double estimatedDegree;

	int isOrphaned;
	int orphanedAt;

	node(int isNodeId){
		this->nodeId = isNodeId;
		this->startTime = 0;
		this->endTime = -1;
		this->estimatedDegree = 0;

		this->isOrphaned = 0;
		this->orphanedAt = -1;
	}

	node(int isNodeId, int isStartTime){
		this->nodeId = isNodeId;
		this->startTime = isStartTime;
		this->endTime = -1;
		this->estimatedDegree = 0;

		this->isOrphaned = 0;
		this->orphanedAt = -1;
	}

	bool addEdge(edge *e){
		if (e->destId == nodeId){
			delete e;
			return false;
		}//no self loops
		bool flag = true;
		for (int i=0;i<adj.size();i++){
			if ((adj[i]->destId == e->destId)){
				if ((adj[i]->endTime == -1)||(adj[i]->endTime > e->startTime)){
					flag = false;
					break;
				}
			}
		}
		if (flag) adj.push_back(e);
		else delete e;  //no multi edges
		return flag;
	}

	/*to be called on destId*/
	edge *findReverseAliveEdge(int srcId){
		for (int i=0;i<adj.size();i++){
			if ((adj[i]->destId == srcId) && (adj[i]->communityId == -1) && (adj[i]->endTime == -1))
				return adj[i];
		}
		return NULL;
	}

	/*to be called on destId*/
	edge *findReverseEpsAliveEdge(int srcId, int commId){
		for (int i=0;i<adj.size();i++){
			if ((adj[i]->destId == srcId) && (adj[i]->communityId == commId) && (adj[i]->endTime == -1))
				return adj[i];
		}
		return NULL;
	}

	void printCommunities(){
		for (int i=0;i<communities.size();i++) cout << ", " << communities[i];
		cout << endl << flush;
	}

	void printEdgeList(){
		cout << "Edge list of " << nodeId << endl;
		for (int i = 0; i < adj.size(); i++)
			cout << "(" << adj[i]->destId  << ", " << adj[i]->startTime << ", " << adj[i]->endTime << ", " << adj[i]->communityId << ") ";
		cout << endl;
	}
	node(){}

	bool hasEdge(int destId){
		bool flag = false;
		for (int i = 0; i < adj.size(); i++){
			if ((adj[i]->destId == destId))
				flag = true;
		}
		return flag;
	}
};

/*A node attaches to a community via this structure.
 It indicates the joining and leaving time of the node in the
 community.*/
struct nodeInCommunity{
	int nodeId;
	int joinTime;
	int leaveTime;

	nodeInCommunity(int isNodeId){
		this->nodeId = isNodeId;
		this->joinTime = -1;
		this->leaveTime = -1;
	}

	nodeInCommunity(int isNodeId, int isJoinTime){
		this->nodeId = isNodeId;
		this->joinTime = isJoinTime;
		this->leaveTime = -1;
	}
};

struct community{
	vector<nodeInCommunity*> nodeList;
	int birthTime;
	int expansionTime;
	int deathTime;
	int contractionTime;
	int nextAvailableTimeSlot;
	bool isAvailable;
	int originFlag; 	//signifies how the community came into being. 0 for birth, 1 for split, 2 for merge
	int deletionFlag;   //signifies how the community stopped being. 0 for death, 1 for split, 2 for merge

	community(){
		this->isAvailable = true;
		this->originFlag = 0;
		this->deletionFlag = 0;
		this->birthTime = 0;
		this->deathTime = -1;
	}

	int indexOfNode(int nodeId){
		for (int i=0;i<nodeList.size();i++)
			if (nodeList[i]->nodeId == nodeId) return i;
		return -1;
	}

	void swapToEnd(int nodeIndex){
		int lastIndex = nodeList.size() - 1;
		nodeInCommunity *nic = nodeList[lastIndex];
		nodeList[lastIndex] = nodeList[nodeIndex];
		nodeList[nodeIndex] = nic;
	}

	void printNodes(){
		for (int i=0;i<nodeList.size();i++) cout << ", " << nodeList[i]->nodeId;
		cout << endl << flush;
	}
};

struct update{
	int updateType; //0: edge delete; 1: node add; 2: edge add; 3: node delete;
	int u;
	int v;
	int t;

	update(int isUpdateType, int isU, int isV, int isT){
		this->updateType = isUpdateType;
		this->u = isU;
		this->v = isV;
		this->t = isT;
	}
};

vector<node*> graph;
vector<community*> communities;
vector<int> nodeMemberships;
BucketSampling *sampler = nullptr;

vector<int> communitySizes;

int numOrphanNodes;
/************PROFILING QUANTITIES************************************/
double averageSplitTime = 0.0;
double averageMergeTime = 0.0;
double averageBirthTime = 0.0;
double averageDeathTime = 0.0;
double averageNodeAddTime = 0.0;
double averageNodeDeleteTime = 0.0;
int nSplits = 0;
int nMerges = 0;
int nBirths = 0;
int nDeaths = 0;
int nAdd = 0;
int nDelete = 0;
/*******************************************************************/

/************************UTILITY FUNCTIONS**************************/
double expectedPowerLaw(int xmin, int xmax, double beta){
	PowerlawDegreeSequence z(xmin,xmax,beta);
	z.run();

	double meanValue = z.getExpectedAverageDegree();
	return meanValue;
}

vector<int> powerLawDegreeSequence(int xmin, int xmax, double beta, int N){
	PowerlawDegreeSequence z(xmin,xmax,beta);
	z.run();

	vector<int> degreeSequence = z.getDegreeSequence(N);
	return degreeSequence;
}

vector<int> powerLawDegreeSequenceSum(int xmin, int xmax, double beta, int sum){
	PowerlawDegreeSequence z(xmin,xmax,beta);
	z.run();

	vector<int> degreeSequence;
	int tempSum = 0;
	while (tempSum < sum){
		int t = z.getDegree();
		if ((tempSum + t) <= sum){
			degreeSequence.push_back(t);
			tempSum += t;
		}
		else{ //whichever is closer
			int d1 = sum - tempSum;
			int d2 = (tempSum + t) - sum;
			if (d2 < d1){
				degreeSequence.push_back(t);
				tempSum += t;
			} else {
				break;
			}
		}
	}
	return degreeSequence;
}

void permuteFY(int *sequence, int size){
	for (int i=0;i<size-1;i++){
		int j = i + rand()%(size-i);
		int temp = sequence[i];
		sequence[i] = sequence[j];
		sequence[j] = temp;
	}
}

int get_next_edge_distance(const double log_cp) {
	return (int) (1 + floor(log(1.0 - (((double) rand())/((double) RAND_MAX))) / log_cp));
}

bool isInCommunityAtT(int t, int communityId, int nodeId){
	community c = *communities[communityId];
	int pos = c.indexOfNode(nodeId);
	bool startFlag = false, endFlag = false;
	if (c.nodeList[pos]->joinTime <= t) startFlag = true;
	if (c.nodeList[pos]->leaveTime == -1) endFlag = true;
	if (c.nodeList[pos]->leaveTime >= t) endFlag = true;
	return (startFlag && endFlag);
}

bool compareUpdate(update i, update j){
	if (i.t < j.t) return true;
	if (i.t > j.t) return false;
	if ((i.u == j.u) && (i.v == j.v)){
		if (i.updateType < j.updateType) return false;
		if (i.updateType > j.updateType) return true;
	}
	if (i.updateType < j.updateType) return true;
	if (i.updateType > j.updateType) return false;

	return false;
}
/*******************************************************************/

/************************STATIC STRUCTURE***************************/
void generateBigraph(){
	nodeMemberships = powerLawDegreeSequence(xmin,xmax,beta1,N1);

	cout << "generateBigraph: CP0" << endl << flush;
	int sumMemberships = 0;
	for (int i = 0; i < N1; i++){
		sumMemberships += nodeMemberships[i];
		//cout << "nM = " << nodeMemberships[i] << endl << flush;
		sampler->addNode(nodeMemberships[i]);
	}

	int newN1 = sampler->getNumberOfNodes();
	for (int i = N1; i < newN1; i++) graph.push_back(new node(i));
	N1 = newN1;

	nodeMemberships.clear();
	nodeMemberships.resize(N1, 0);

	cout << "sampler has added " << sampler->getNumberOfNodes() << endl << flush;

	/*Generate power law degree sequences*/
	communitySizes = powerLawDegreeSequenceSum(mmin,mmax,beta2,sumMemberships);
	N2 = communitySizes.size();
	int actualSumOfMemberships = 0;

	for (int i=0;i<N2;i++){
		actualSumOfMemberships += communitySizes[i];
		communities.push_back(new community());
		if (actualSumOfMemberships > sampler->getSumOfDesiredMemberships()){
			break;
		}
	}
	N2 = communities.size();

	cout << "GenerateBigraph : CP1" << endl << flush;
	for (int i=0;i<N2;i++){
		communities[i]->birthTime = 0;
		int wanted_size = communitySizes[i];
		std::vector<int> communityNodes = sampler->birthCommunityNodes(wanted_size);
		cout << "Generating community : " << i << ", " << communityNodes.size() << ", " << wanted_size << endl << flush;
		for (int u : communityNodes) {
			sampler->assignCommunity(u);
			++nodeMemberships[u];
			++currMem;
			graph[u]->communities.push_back(i);
			(communities[i]->nodeList).push_back(new nodeInCommunity(u,0));
		}
		initialMem += communityNodes.size();
	}
}

void generateEdgesForCommunity(int commIndex){
	vector<int> nodesInCommunity;
	for (int i=0;i<N1;i++){
		for (int j=0;j<graph[i]->communities.size();j++)
			if (graph[i]->communities[j] == commIndex){
				nodesInCommunity.push_back(i);
				break;
			}
	}
	int numberOfNodes = nodesInCommunity.size();
	if (numberOfNodes <= 1) return;

	double probNew = alpha/pow(numberOfNodes, gamma_1);

	if (numberOfNodes == 2) probNew = 1;
	if (probNew > 1.0) probNew = 1.0;
	/*Copied from Networkit implementation of Batagelj Brandes*/
	const double log_cp = log(1.0 - probNew); // log of counter probability
	// create edges
	int curr = 1;
	int next = -1;
	while (curr < numberOfNodes) {
		// compute new step length
		next += get_next_edge_distance(log_cp);
		// check if at end of row
		while ((next >= curr) && (curr < numberOfNodes)) {
			// adapt to next row
			next = next - curr;
			curr++;
		}
		// insert edge
		if (curr < numberOfNodes) {
			int a = nodesInCommunity[curr];
			int b = nodesInCommunity[next];

			edge *fwd = new edge(a,b,commIndex);
			fwd->startTime = 0;
			bool flagR = graph[a]->addEdge(fwd);
			if (flagR){
				edge *bwd = new edge(b,a,-1);    //communityId = -1 for a reverse edge
				bwd->startTime = 0;
				flagR = graph[b]->addEdge(bwd);
				if (!flagR){
					cout << "1.ERROR HERE!" << endl << flush;
					exit(0);
				}
			}

		}
	}

	double expectedDegree = (numberOfNodes-1)*probNew;
	for (int i=0;i<nodesInCommunity.size();i++)
		graph[i]->estimatedDegree += ((int)(ceil(expectedDegree)));

}

void generateEpsCommunity(){
	double epsilon = eps/N1;
	int numEdges = (int) floor(epsilon * N1 * (N1-1) * 0.5);
	for (int i=0;i<numEdges;i++){
		int sourceNode1 = rand()%N1;
		int destNode1 = rand()%N1;
		int sourceNode2 = rand()%N1;
		int destNode2 = rand()%N1;

		int switchTime = T/2 + (rand()%200 - 100);
		if (switchTime < 0) switchTime = 0;
		if (sourceNode1 != destNode1){
			edge *fwd = new edge(sourceNode1,destNode1,-2);    //communityId = -2 for an external edge
			fwd->startTime = 0;
			fwd->endTime = switchTime;
			bool flagR = graph[sourceNode1]->addEdge(fwd);
			if (flagR){
				edge *bwd = new edge(destNode1,sourceNode1,-4);
				bwd->startTime = 0;
				bwd->endTime = switchTime;
				flagR = graph[destNode1]->addEdge(bwd);
				if (!flagR){
					cout << "2.ERROR HERE!" << endl << flush;
					exit(0);
				}
			}
		}
		if (sourceNode2 != destNode2){
			edge *fwd = new edge(sourceNode2,destNode2,-2);    //communityId = -2 for an external edge
			fwd->startTime = switchTime + 1;
			fwd->endTime = -1;
			bool flagR = graph[sourceNode2]->addEdge(fwd);
			if (flagR){
				edge *bwd = new edge(destNode2,sourceNode2,-4);
				bwd->startTime = switchTime + 1;
				bwd->endTime = -1;
				flagR = graph[destNode2]->addEdge(bwd);
				if (!flagR){
					cout << "3.ERROR HERE!" << endl << flush;
					exit(0);
				}
			}
		}
	}
}

bool isSanity = true;
void sanityCheck(){
	for (int i=0;i<N1;i++){
		for (int j=0; j< graph[i]->communities.size(); j++){
			int cIndex = graph[i]->communities[j];
			community *c = communities[cIndex];
			if (c->indexOfNode(i) == -1) isSanity = false;
		}
	}
}

bool sanityCheck2(){
	for (int i=0;i<N1;i++){
		for (int j = 0; j< graph[i]->adj.size(); j++){
			edge *e = graph[i]->adj[j];
			if ((e->startTime > e->endTime) && (e->startTime != -1) && (e->endTime != -1))
				return false;
		}
	}
	return true;
}

bool sanityCheck3(){
	for (int i=0;i<N1;i++){
		for (int j = 0; j < graph[i]->adj.size(); j++){
			int dj = graph[i]->adj[j]->destId, sj = graph[i]->adj[j]->startTime, ej = graph[i]->adj[j]->endTime;
			for (int k = j+1; k < graph[i]->adj.size(); k++){
				int dk = graph[i]->adj[k]->destId, sk = graph[i]->adj[k]->startTime, ek = graph[i]->adj[k]->endTime;
				if (dj == dk){
					//intervals should be disjoint
					bool flag = false;
					if ((sj >= ek) && (ek != -1)) flag = true;
					if ((sk >= ej) && (ej != -1)) flag = true;
					if (!flag){
						cout << "( "<< sj << ", " << ej << " ), ( " << sk << ", " << ek << " )" << ", " << graph[i]->adj[j]->communityId << ", " << graph[i]->adj[k]->communityId << endl;
						return false;
					}
				}
			}
		}
	}
	return true;
}

void generateStaticStructure(){
	for (int i=0;i<N1;i++) graph.push_back(new node(i));
	cout << "Generating bigraph... " << endl << flush;
	generateBigraph();
	cout << "Generated bigraph... " << endl << flush;
	int numCommunities = communitySizes.size();
	cout << "numCommunities = " << numCommunities << endl << flush;
	for (int i=0;i<numCommunities;i++){
		generateEdgesForCommunity(i);
		if (i % 10 == 0) cout << "done " << i << " comms " << endl << flush;
	}
	generateEpsCommunity();
	sanityCheck();
}
/*******************************************************************/

/*************************COMMUNITY EVENTS***************************/
void deathCommunity(int commIndex, int timeslot){
	community *c = communities[commIndex];
	c->isAvailable = false;
	c->nextAvailableTimeSlot = -1;
	c->deathTime = timeslot;
	c->contractionTime = timeslot - rand()%10 + 5; //contraction time between 5 and 14 s
	if (c->contractionTime < 0) c->contractionTime = 0;
	int P = c->nodeList.size();
	int coreSize = (int)(0.1*P);
	if (coreSize < mmin) coreSize = mmin;
	for (int i=0 ; i < c->nodeList.size() ; i++){
        if (c->nodeList[i]->leaveTime > -1) continue;
		if (i < coreSize) c->nodeList[i]->leaveTime = c->deathTime;
		else c->nodeList[i]->leaveTime = (int) floor(c->deathTime - (((double) (i - coreSize))/(P - coreSize))*(c->deathTime - c->contractionTime));
	}

	for (int i=0;i<c->nodeList.size();i++){
		int u = c->nodeList[i]->nodeId;
		if (graph[u]->endTime > -1) continue;   //this vertex is already dead, and so has left the community
        	nodeMemberships[u] -= 1; currMem -= 1;
		sampler->leaveCommunity(u);
        	if (nodeMemberships[u] == 0){
            		numOrphanNodes += 1;
            		cout << "Orphaned " << u << endl;
	    		graph[u]->isOrphaned = 1;
	    		graph[u]->orphanedAt = timeslot;
			//delete all epsilon edges.
			for (int j = 0; j < graph[u]->adj.size(); j++){
				if (graph[u]->adj[j]->endTime == -1 && graph[u]->adj[j]->communityId == -2){
					int dj = graph[u]->adj[j]->destId;
					edge *reverseEdge = graph[dj]->findReverseEpsAliveEdge(u, -4);
					if ((reverseEdge) && (reverseEdge->endTime == -1) && (reverseEdge->startTime <= timeslot)){
						reverseEdge->endTime = timeslot;
						graph[u]->adj[j]->endTime = timeslot;
					}
				}
				else if (graph[u]->adj[j]->endTime == -1 && graph[u]->adj[j]->communityId == -4){
					int dj = graph[u]->adj[j]->destId;
					edge *reverseEdge = graph[dj]->findReverseEpsAliveEdge(u, -2);
					if ((reverseEdge) && (reverseEdge->endTime == -1) && (reverseEdge->startTime <= timeslot)){
						reverseEdge->endTime = timeslot;
						graph[u]->adj[j]->endTime = timeslot;
					}
				}
			}

		}

		for (int j=0; j< graph[u]->adj.size();j++){
			if (graph[u]->adj[j]->communityId == commIndex && graph[u]->adj[j]->endTime == -1){
				//generate end time with exponential
				int p_u = i;
				int p_v = c->indexOfNode(graph[u]->adj[j]->destId);
				int uLeaveTime = c->nodeList[p_u]->leaveTime;
				int vLeaveTime = c->nodeList[p_v]->leaveTime;
				int minLeaveTime = uLeaveTime;
				if (minLeaveTime > vLeaveTime) minLeaveTime = vLeaveTime;
				int endLimit = minLeaveTime;
				graph[u]->adj[j]->generateEndExp(endLimit);
				//set end time of reverse edge here
				int dj = graph[u]->adj[j]->destId;
				edge *reverseEdge = graph[dj]->findReverseAliveEdge(u);
				if ((reverseEdge)&&(reverseEdge->endTime == -1)&&(reverseEdge->startTime <= graph[u]->adj[j]->endTime)){
		                    reverseEdge->endTime = graph[u]->adj[j]->endTime;
				}
			}
		}
	}
}

void birthCommunity(int timeslot){
    	int includedNode = -1;
	//draw community size
	PowerlawDegreeSequence z(mmin,mmax,beta2);
	z.run();
	int commSize = z.getDegree();
    	if (commSize < mmin) return;

	std::vector<int> commNodes = sampler->birthCommunityNodes(commSize);
	int commIndex = communities.size();
	//cout << "Birthing community " << commIndex << " at " << timeslot << endl;
	//cout << "N2 = " << N2 << endl;
	communities.push_back(new community());
	N2 += 1;
	community *c = communities[commIndex];
	c->birthTime = timeslot;
	c->expansionTime = timeslot + rand()%10 + 5; //contraction time between 5 and 14 s

	for (int i=0;i<commSize;i++){
		sampler->assignCommunity(commNodes[i]);
		c->nodeList.push_back(new nodeInCommunity(commNodes[i]));
		graph[commNodes[i]]->communities.push_back(commIndex);
		nodeMemberships[commNodes[i]] += 1; currMem += 1;
		if (nodeMemberships[commNodes[i]] == 1){
			int u = commNodes[i];
			numOrphanNodes -= 1;
			graph[u]->isOrphaned = 0;
			//recover all epsilon edges
			for (int j = 0; j < graph[u]->adj.size(); j++){
				int dj = graph[u]->adj[j]->destId;
				if (graph[u]->adj[j]->endTime == graph[u]->orphanedAt
					&& graph[u]->adj[j]->communityId == -2
					&& graph[dj]->isOrphaned == 0){
					edge *fwd = new edge(u, dj, -2);
					edge *bwd = new edge(dj, u, -4);
					fwd->startTime = bwd->startTime = timeslot;
					if (graph[u]->addEdge(fwd)) graph[dj]->addEdge(bwd);
				}
				else if (graph[u]->adj[j]->endTime == graph[u]->orphanedAt
					&& graph[u]->adj[j]->communityId == -4
					&& graph[dj]->isOrphaned == 0){
					edge *fwd = new edge(u, dj, -4);
					edge *bwd = new edge(dj, u, -2);
					fwd->startTime = bwd->startTime = timeslot;
					if (graph[u]->addEdge(fwd)) graph[dj]->addEdge(bwd);
				}
			}
			graph[u]->orphanedAt = -1;
		}
	}


	//generate internal edges
	//double probNew = 2*prob/(numberOfNodes-1);
	double probNew = alpha/pow(commSize,gamma_1);
	if ((probNew > 1.0) || (commSize == 2)) probNew = 1;
	int coreSize = (int)(0.1*commSize);
	if (coreSize < mmin) coreSize = mmin;
	for (int i=0; i<commSize; i++){
		if (i < coreSize) c->nodeList[i]->joinTime = c->birthTime;
		else c->nodeList[i]->joinTime = (int) (floor(c->birthTime + (((double) (i - coreSize))/(commSize - coreSize))*(c->expansionTime - c->birthTime)));
	}

	/*Copied from Networkit implementation of Batagelj Brandes*/
	const double log_cp = log(1.0 - probNew); // log of counter probability
	// create edges
	int curr = 1;
	int next = -1;
	c->nextAvailableTimeSlot = c->expansionTime;
	cout << "Birth community: " << c->expansionTime << endl;

	while (curr < commSize) {
		// compute new step length
		next += get_next_edge_distance(log_cp);
		// check if at end of row
		while ((next >= curr) && (curr < commSize)) {
			// adapt to next row
			next = next - curr;
			curr++;
		}
		// insert edge
		if (curr < commSize) {
			int a = commNodes[curr];
			int b = commNodes[next];
			edge *fwd = new edge(a,b,commIndex);
			bool flag = graph[a]->addEdge(fwd);
			if (!flag) continue;
			edge *bwd = new edge(b,a,-1);    //communityId = -1 for a reverse edge
			bwd->startTime = timeslot;
			flag = graph[b]->addEdge(bwd);

			//generate start time
			int p_u = c->indexOfNode(a);
			int p_v = c->indexOfNode(b);
			int uJoinTime = c->nodeList[p_u]->joinTime;
			int vJoinTime = c->nodeList[p_v]->joinTime;
			int maxJoinTime = uJoinTime;
			if (maxJoinTime < vJoinTime) maxJoinTime = vJoinTime;
			int startLimit = maxJoinTime;
			fwd->generateStartExp(startLimit);
			bwd->startTime = fwd->startTime;
		}
	}
	//mark busy
	c->isAvailable = false;
}

void generateEdgesForSplitCommunity(int commIndex, int timeslot){
	community *c = communities[commIndex];
	int numberOfNodes = c->nodeList.size();
	if (numberOfNodes <= 1) return;

	double p1 = alpha/pow(numberOfNodes, gamma_1);
	double p0 = 0;
	for (int i=0; i<numberOfNodes; i++){
		int u = c->nodeList[i]->nodeId;
		for (int j=0; j < graph[u]->adj.size(); j++){
			if (graph[u]->adj[j]->communityId == commIndex)
				p0 += 1;
		}
	}
	p0 /= (numberOfNodes*(numberOfNodes-1))/2.0;
	if (p1 <= p0) return;
	double probNew = (p1-p0)*(1+p0);
	if (numberOfNodes == 2) probNew = 1;
	if (probNew > 1.0) probNew = 1.0;
	/*Copied from Networkit implementation of Batagelj Brandes*/
	const double log_cp = log(1.0 - probNew); // log of counter probability
	// create edges
	int curr = 1;
	int next = -1;
	while (curr < numberOfNodes) {
		// compute new step length
		next += get_next_edge_distance(log_cp);
		// check if at end of row
		while ((next >= curr) && (curr < numberOfNodes)) {
			// adapt to next row
			next = next - curr;
			curr++;
		}
		// insert edge
		if (curr < numberOfNodes) {
			int a = c->nodeList[curr]->nodeId;
			int b = c->nodeList[next]->nodeId;
			edge *fwd = new edge(a,b,commIndex);
			fwd->startTime = timeslot;
			bool flag = graph[a]->addEdge(fwd);
			if (flag){
				edge *bwd = new edge(b,a,-1);    //communityId = -1 for a reverse edge
				bwd->startTime = timeslot;
				flag = graph[b]->addEdge(bwd);
				if (!flag){
					cout << "5.ERROR HERE!" << endl << flush;
					graph[a]->printEdgeList();
					graph[b]->printEdgeList();
					exit(0);
				}
			}
		}
	}
}

void splitCommunity(int commIndex, int timeslot){

	//cout << "N2 = " << N2 << endl << flush;
	community *c = communities[commIndex];

	c->isAvailable = false;
	c->nextAvailableTimeSlot = -1;  //community dead, will never be available again
	c->deletionFlag = 1;
	double splitPoint = minSplitRatio + (((double) rand())/((double) RAND_MAX))*(1.0 - 2.0*minSplitRatio);
	int splitBorder = (floor) (splitPoint * c->nodeList.size());
	int L = splitBorder+1;
	int R = c->nodeList.size() - L;


	//create the two smaller communities
	int c1Index = communities.size();
	communities.push_back(new community());
	community *c1 = communities[c1Index];
	c1->originFlag = 1;
	c1->isAvailable = false;


	int c2Index = communities.size();
	communities.push_back(new community());
	community *c2 = communities[c2Index];
	c2->originFlag = 1;
	c2->isAvailable = false;
	cout << "Splitting community " << commIndex << ", of size = " << c->nodeList.size() << ", into " << c1Index << " and " << c2Index << " at " << timeslot << ", originFlag = " << c->originFlag << endl << flush;

	int *leaveTimes = new int[c->nodeList.size()];
	int numNodesInvolved = 0;
	for (int i=0;i < (c->nodeList.size());i++){
		int u = (c->nodeList[i])->nodeId;
		if (graph[u]->endTime > -1){
			leaveTimes[i] = -1;
			continue;
		}
		numNodesInvolved += 1;
		c->nodeList[i]->leaveTime = timeslot;
		nodeInCommunity *nic = new nodeInCommunity(u,timeslot);
		if (i>splitBorder){
			leaveTimes[i] = (int) floor(timeslot + (((double) (R+L-i))/R)*timeToSplit);
			c2->nodeList.push_back(nic);
			graph[u]->communities.push_back(c2Index);
		}
		else{
			leaveTimes[i] = (int) floor(timeslot + (((double) i)/L)*timeToSplit);
			c1->nodeList.push_back(nic);
			graph[u]->communities.push_back(c1Index);
		}
	}
	//generate edge end times
	int maxEndTime = 0;

	for (int i=0;i<=splitBorder;i++){	//for all nodes in c1
		int u = c->nodeList[i]->nodeId;
		for (int j=0;j<graph[u]->adj.size();j++){
			if (graph[u]->adj[j]->communityId != commIndex || graph[u]->adj[j]->endTime != -1) continue;    //is a different edge, no need to generate end time

			int v = graph[u]->adj[j]->destId;
			if (graph[v]->endTime > -1) continue;
			int indexOfv = c->indexOfNode(v);
			if (indexOfv > splitBorder){	//if v is in c2
				int endLimit = (int) floor(0.5*leaveTimes[i] + 0.5*leaveTimes[indexOfv]);
				if (endLimit < timeslot){
					cout << "1.Now this problem is there!" << endl;
				}
				graph[u]->adj[j]->generateEndExp(endLimit);
				//set endTime of reverse edge here
				int dj = graph[u]->adj[j]->destId;
				edge *reverseEdge = graph[dj]->findReverseAliveEdge(u);
				if (reverseEdge){
					if (reverseEdge->endTime == -1)
						reverseEdge->endTime = graph[u]->adj[j]->endTime;
				}
				if (maxEndTime < graph[u]->adj[j]->endTime) maxEndTime = graph[u]->adj[j]->endTime;
			}
			else graph[u]->adj[j]->communityId = c1Index;
		}
	}
	/*This loop is required because we are checking if the edge processed has community id = commIndex, so reverse edges
	will not be processed in the previous one.*/
	for (int i=splitBorder+1;i<numNodesInvolved;i++){	//for all nodes in c2
		int u = c->nodeList[i]->nodeId;
		for (int j=0;j<graph[u]->adj.size();j++){
			if (graph[u]->adj[j]->communityId != commIndex || graph[u]->adj[j]->endTime != -1) continue;    //is a different edge, no need to generate end time
			int v = graph[u]->adj[j]->destId;
			if (graph[v]->endTime > -1) continue;   //v is already dead
			int indexOfv = c->indexOfNode(v);
			if (indexOfv <= splitBorder && indexOfv >=0){ //if v is in c1
				int endLimit = (int) floor(0.5*leaveTimes[indexOfv] + 0.5*leaveTimes[i]);
				if (endLimit < timeslot){
					cout << "2. Now this problem is there!" << endl;
				}
				graph[u]->adj[j]->generateEndExp(endLimit);
				//set end time of reverse edge here
				int dj = graph[u]->adj[j]->destId;
				edge *reverseEdge = graph[dj]->findReverseAliveEdge(u);
				if ((reverseEdge)&&(reverseEdge->endTime==-1))
					reverseEdge->endTime = graph[u]->adj[j]->endTime;

				if (maxEndTime < graph[u]->adj[j]->endTime) maxEndTime = graph[u]->adj[j]->endTime;
			}
			else graph[u]->adj[j]->communityId = c2Index;
		}
	}

	c1->birthTime = timeslot + 1;
	c2->birthTime = timeslot + 1;
	c->deathTime = timeslot;
    cout << "Split community: " << maxEndTime + 1 << endl;
	c1->nextAvailableTimeSlot = maxEndTime+1;
	c2->nextAvailableTimeSlot = maxEndTime+1;
	generateEdgesForSplitCommunity(c1Index,timeslot);
	generateEdgesForSplitCommunity(c2Index,timeslot);
	N2 += 2;
	delete [] leaveTimes;
}

void mergeCommunities(int c1Index, int c2Index, int timeslot){

	community *c1 = communities[c1Index];
	community *c2 = communities[c2Index];


	c1->isAvailable = false;
	c1->nextAvailableTimeSlot = -1; //community dead, will never be available again
	c1->deletionFlag = 2;

	c2->isAvailable = false;
	c2->nextAvailableTimeSlot = -1;
	c2->deletionFlag = 2;

	int cIndex = communities.size();
	communities.push_back(new community());
	community *c = communities[cIndex];
	cout << "Merging communities " << c1Index << " and " << c2Index << " into " << cIndex << " at " << timeslot << endl;
	//cout << "N2 = " << N2 << endl << flush;
	c->originFlag = 2;
	c->isAvailable = false;
	int L = c1->nodeList.size();
	int R = c2->nodeList.size();

	for (int i=0;i<c1->nodeList.size();i++){
		int u = c1->nodeList[i]->nodeId;
		if (graph[u]->endTime > -1) continue;
		c1->nodeList[i]->leaveTime = timeslot + timeToMerge;
		nodeInCommunity *nic = new nodeInCommunity(u, timeslot + timeToMerge + 1);
		c->nodeList.push_back(nic);
		graph[u]->communities.push_back(cIndex);

		if (c2->indexOfNode(u) != -1){
			nodeMemberships[u] -= 1;
			sampler->leaveCommunity(u);
			currMem -= 1;
			continue;
		}

		nodeInCommunity *nic1 = new nodeInCommunity(u, (int) floor(timeslot + (1-(((double) i)/c1->nodeList.size()))*timeToMerge));
		nic1->leaveTime = timeslot + timeToMerge;
		c2->nodeList.push_back(nic1);
		graph[u]->communities.push_back(c2Index);
		//cout << "Inserted " << u << " to c2" << endl;
	}
	for (int i=0;i<R;i++){
		int u = c2->nodeList[i]->nodeId;
		if (graph[u]->endTime > -1) continue;
		c2->nodeList[i]->leaveTime = timeslot + timeToMerge;
		int posU = c->indexOfNode(u);
		if (posU != -1) continue;
		nodeInCommunity *nic = new nodeInCommunity(u, timeslot + timeToMerge + 1);
		c->nodeList.push_back(nic);
		graph[u]->communities.push_back(cIndex);

		if (c1->indexOfNode(u) != -1) continue;

		nodeInCommunity *nic1 = new nodeInCommunity(u, (int) floor(timeslot + (((double) i)/c2->nodeList.size())*timeToMerge));
		nic1->leaveTime = timeslot + timeToMerge;
		c1->nodeList.push_back(nic1);
		graph[u]->communities.push_back(c1Index);
		//cout << "Inserted " << u << " to c1" << endl;
	}

	/*number of edges to be inserted between nodes of erstwhile c1 and c2
	* depends on the balance number of edges that remains from
	* the required number in c minus the existing number in c1 + c2*/

	int maxStartTime=0;
	int numberOfNodes = c->nodeList.size();
	double prob = alpha/pow(numberOfNodes, gamma_1);

	double m1 = 0, m2 = 0, m = prob*(numberOfNodes)*(numberOfNodes-1)*0.5;
	for (int i=0; i< c1->nodeList.size();i++){
        int c1u = c1->nodeList[i]->nodeId;
		for (int j=0;j<graph[c1u]->adj.size();j++)
			if (graph[c1u]->adj[j]->communityId == c1Index){
				m1 += 1;
				graph[c1u]->adj[j]->communityId = cIndex;
			}
    }
	for (int i=0; i< c2->nodeList.size();i++){
        int c2u = c2->nodeList[i]->nodeId;
		for (int j=0;j<graph[c2u]->adj.size();j++)
			if (graph[c2u]->adj[j]->communityId == c2Index){
				m2 += 1;
				graph[c2u]->adj[j]->communityId = cIndex;
			}
    }
	double balanceEdges = m - m1 - m2;
	double probNew = balanceEdges/(L*R);
	//generate edges and their start times

	if ((probNew > 1.0) || (numberOfNodes == 2)) probNew = 1;
	if (probNew < 0.0) probNew = 0;
	//cout << "L = " << L << ", R = " << R << ", timeToMerge = " << timeToMerge << endl << flush;
	for (int i=0;i<L;i++){
		int u = c1->nodeList[i]->nodeId;
		if (graph[u]->endTime > -1) continue;
		int posU = c2->indexOfNode(u);
		for (int j=0;j<R;j++){
			int v = c2->nodeList[j]->nodeId;
			if (graph[v]->endTime > -1) continue;
			/*bool flag = true;
			for (int l=0;l<graph[u]->adj.size();l++)
				if (graph[u]->adj[l]->destId == v){
					flag = false;
					break;
				}
			if (!flag) continue;*/
			double r = ((double) rand())/((double) RAND_MAX);
			if (r<=probNew){
				edge *fwd = new edge(u,v,cIndex);
				int posV = c1->indexOfNode(v);

				int startLimit = (int) floor(0.5*c1->nodeList[posV]->joinTime + 0.5*c2->nodeList[posU]->joinTime);
				if (startLimit < timeslot){	//this can happen because of a pre-existing overlap between c1 and c2
					//cout << "i = " << i << ", j = " << j << ", posU = " << posU << ", posV = " << posV << ", startLimit = " << startLimit << ", c1->nodeList[posV].joinTime = " << c1->nodeList[posV].joinTime << ", c2->nodeList[posU].joinTime = " << c2->nodeList[posU].joinTime <<  endl;
					startLimit = timeslot;
				}
				fwd->generateStartExp(startLimit);
				if (maxStartTime < fwd->startTime) maxStartTime = fwd->startTime;
				bool flag = graph[u]->addEdge(fwd);
				if (flag){
					edge *bwd = new edge(v,u,-1);
					bwd->startTime = fwd->startTime;
					flag = graph[v]->addEdge(bwd);
					if (!flag){
						cout << "6.ERROR HERE!" << endl << flush;
						exit(0);
					}
				}
			}
		}
	}
	if (maxStartTime < timeslot) cout << "MaxStartTime Issue: probNew = " << probNew << endl;
    if (maxStartTime < timeslot) maxStartTime = timeslot;
	c->birthTime = timeslot + timeToMerge + 1;
	c1->deathTime = timeslot + timeToMerge;
	c2->deathTime = timeslot + timeToMerge;
    cout << "Merge community: " << maxStartTime + 1 << endl;
	c->nextAvailableTimeSlot = maxStartTime+1;
	N2 += 1;


	/*cout << "End: c2s nodeList" << endl;
	for (int i=0;i<R;i++)
		cout << "id: " << c2->nodeList[i].nodeId << ", join = " << c2->nodeList[i].joinTime << ", leave = " << c2->nodeList[i].leaveTime << endl;

	cout << "Finished Merging" << endl;
	*/

}
/*********************************************************************/

/**************************VERTEX EVENTS****************************/
void addNode(int timeslot){
	PowerlawDegreeSequence z(xmin,xmax,beta1); //get community membership
	z.run();
	int nodeMembership = z.getDegree();
	sampler->addNode(nodeMembership);
	N1 = sampler->getNumberOfNodes();
	//nodeMemberships.push_back((int)(floor(1.2*nodeMembership)));
	while (graph.size() < N1) {
	    const int u = graph.size();
	    graph.push_back(new node(u,timeslot));
	    nodeMemberships.push_back(0);
	    for (int i=0;i<eps;i++){
		    int v = rand()%(N1-1);
		    while (graph[v]->endTime != -1)
			    v = rand()%(N1-1);

		    edge *fwd = new edge(u,v,-2);
		    fwd->startTime = timeslot;
		    bool flag = graph[u]->addEdge(fwd);
		    if (flag){
			    edge *bwd = new edge(v,u,-1);    //communityId = -1 for a reverse edge
			    bwd->startTime = timeslot;
			    flag = graph[v]->addEdge(bwd);
			    if (!flag){
				    cout << "7.ERROR HERE!" << endl << flush;
				    exit(0);
			    }
		    }
	    }
	}
}

void deleteNode(int nodeIndex, int timeslot){
	sampler->removeNode(nodeIndex);

	nodeMemberships[nodeIndex] = 0; //no future births will involve this node
	for (int i=0;i<graph[nodeIndex]->communities.size();i++){
		int communityId = graph[nodeIndex]->communities[i];
		community c = *communities[communityId];
		int pos = c.indexOfNode(nodeIndex);
		if (c.nodeList[pos]->leaveTime == -1)
			c.nodeList[pos]->leaveTime = timeslot;
		int numNodesAlive = c.nodeList.size();
		for (int j=0;j<c.nodeList.size();j++)
			if (c.nodeList[j]->leaveTime > -1) numNodesAlive -= 1;
		if (numNodesAlive < mmin){
			for (int j=0;j<c.nodeList.size();j++)
				c.nodeList[j]->leaveTime = timeslot;
			c.isAvailable = false;
			c.nextAvailableTimeSlot = -1;
			c.deathTime = timeslot;
		}
	}
	for (int i=0;i<graph[nodeIndex]->adj.size();i++){
		if (graph[nodeIndex]->adj[i]->endTime == -1) graph[nodeIndex]->adj[i]->endTime = timeslot;
	}
	graph[nodeIndex]->endTime = timeslot;
}

/*******************************************************************/
void generateEvent(int timeslot){
	struct timespec start,finish;
	int z;
	z = rand()%50 + 1;
	if (z==0){
		z = rand()%4;
		if (z>0){
			//insert node
			//cout << "Starting event add node" << endl;
			clock_gettime(CLOCK_MONOTONIC,&start);
			addNode(timeslot);
			clock_gettime(CLOCK_MONOTONIC,&finish);
			averageNodeAddTime += (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
			nAdd += 1;
			//cout << "Finished event add node" << endl;
		}
		else{
			clock_gettime(CLOCK_MONOTONIC, &start);
			int u = rand()%N1;
			while (graph[u]->endTime > -1) u = rand()%N1;
			deleteNode(u,timeslot);
			clock_gettime(CLOCK_MONOTONIC,&finish);
			averageNodeDeleteTime += (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
			nDelete += 1;
		}
		return;
	}

	z = rand()%4;
	double zz = (1.0 * rand())/RAND_MAX ;

	if (zz <= probD){
		//draw an available community at random
		clock_gettime(CLOCK_MONOTONIC,&start);
		int commToDelete = -1, numCandidates = 0;
		for (int i=0;i<N2;i++){
			if (communities[i]->isAvailable){
				numCandidates += 1;
				if (rand()%numCandidates == 0)
					commToDelete = i;
			}
		}

		if (commToDelete != -1) deathCommunity(commToDelete,timeslot);
		cout << "Deleting Community " << commToDelete << " at " << timeslot << endl;
		clock_gettime(CLOCK_MONOTONIC,&finish);
		averageDeathTime += (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
		nDeaths += 1;
	}
	else if (zz <= (probD + probB)){
		clock_gettime(CLOCK_MONOTONIC,&start);
		cout << "Birthing Community " << N2 << " at " << timeslot << endl;
		birthCommunity(timeslot);
		cout << "Community Born " << N2 << " at " << timeslot << endl << flush;
		clock_gettime(CLOCK_MONOTONIC,&finish);
		for (int i = 0; i < N2; i++) communities[i]->isAvailable;
		cout << "Yo. Community Born" << endl << flush;
		averageBirthTime += (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
		nBirths += 1;
	}
	else if (zz <= (probD + probB + probS)){
		int minSplitSize = (int)(ceil(mmin/minSplitRatio));
		clock_gettime(CLOCK_MONOTONIC,&start);
		int commToSplit = -1, numCandidates = 0;
		for (int i=0;i<N2;i++){
			if (communities[i]->isAvailable && (communities[i]->nodeList).size() > (minSplitSize)){
				numCandidates += 1;
				if (rand()%numCandidates == 0)
					commToSplit = i;
			}
		}
		if (commToSplit !=-1){
			splitCommunity(commToSplit,timeslot);
		}
		else{
		int numAvailableCommunities = 0;
            for (int i =0; i< N2; i++) if (communities[i]->isAvailable) numAvailableCommunities++;
            cout << "Attempted split, Could not find suitable community in " << numAvailableCommunities << " to split at " << timeslot << endl;
		}
		clock_gettime(CLOCK_MONOTONIC,&finish);
		averageSplitTime += (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
		nSplits += 1;
	}
	else{ //merge
		clock_gettime(CLOCK_MONOTONIC,&start);
		int comm[2], numCandidates = 0;
		for (int i=0;i<N2;i++){
			if (!communities[i]->isAvailable) continue;
			if (((communities[i]->nodeList).size()) < mmax/2){
				numCandidates += 1;
				int r = rand()%numCandidates;
				if (r==0) comm[0] = i;
			}
		}
		numCandidates = 0;
		int commSize = (communities[comm[0]]->nodeList).size();
		for (int j=0; j<N2;j++){
			if (comm[0]==j) continue;
			if (!communities[j]->isAvailable) continue;
			if (((communities[j]->nodeList).size() + commSize)<=mmax){
				numCandidates += 1;
				int r = rand()%numCandidates;
				if (r==0){
					comm[1] = j;
				}
			}
		}
		if (numCandidates >= 2){
			mergeCommunities(comm[0],comm[1],timeslot);
		}
		else{
            int numAvailableCommunities = 0;
            for (int i =0; i< N2; i++) if (communities[i]->isAvailable) numAvailableCommunities++;
            cout << "Attempted merge, but no suitable candidates of " << numAvailableCommunities << " at " << timeslot << endl;
		}
		clock_gettime(CLOCK_MONOTONIC,&finish);
		averageMergeTime += (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
		nMerges += 1;
	}
	cout << "Finished event " << endl << flush;
}

void randomlyPerturb(community* c, int timeslot, int commIndex){
	//sample a node
	int numberOfNodes = c->nodeList.size();
	int nodeIndex = rand() % numberOfNodes;
	int nodeId = c->nodeList[nodeIndex]->nodeId;
	//sample an edge of this node - to be deleted
	edge *eToD = NULL;
	int numEdges = 0;
	for (int i = 0; i < graph[nodeId]->adj.size(); i++){
		edge *e = graph[nodeId]->adj[i];
		if ((e->endTime == -1)&&(e->communityId == commIndex)&&(e->startTime < timeslot)){
			numEdges += 1;
			int r = rand()%numEdges;
			if (r==0) eToD = e;
		}
	}
	if (eToD){
		if ((eToD->startTime > eToD->endTime) && (eToD->startTime != -1) && (eToD->endTime != -1)){
			cout << "communityId = " << eToD->communityId << ", startTime = " << eToD->startTime << ", endTime = " << eToD->endTime << endl;
			exit(0);
		}
	}
	else return;

	//sample a node pair to be inserted as edge
	int sNode = eToD->sourceId;
	int dNode = -1;
	int numCandidateNodes = 0;
	for (int i = 1; i < numberOfNodes; i++){
		int d = c->nodeList[i]->nodeId;
		if (d == sNode) continue;
		if (graph[sNode]->hasEdge(d) == false){
			numCandidateNodes += 1;
			int r = rand()%(numCandidateNodes);
			if (r==0) dNode = d;
		}
	}
	if (dNode == -1) return;

	//cout << "Inserting Edge " << commIndex << ":" << sNode << ", " << dNode << ", " << timeslot << endl << flush;
	edge *fwd = new edge(sNode,dNode,commIndex);
	fwd->startTime = timeslot;
	bool flagR = graph[sNode]->addEdge(fwd);
	if (flagR){
		edge *bwd = new edge(dNode, sNode, -1);
		bwd->startTime = timeslot;
		flagR = graph[dNode]->addEdge(bwd);
		if (!flagR){
			cout << "8.ERROR HERE!" << endl << flush;
			graph[sNode]->printEdgeList();
			graph[dNode]->printEdgeList();
			exit(0);
		}
	}
	eToD->endTime = timeslot;
	edge *reverseEdge = graph[eToD->destId]->findReverseAliveEdge(eToD->sourceId);
	if ((reverseEdge)&&(reverseEdge->endTime == -1)) reverseEdge->endTime = timeslot;
}

void generateTimeslot(int timeslot){
	double z = ((double) rand())/RAND_MAX;  //10 % probability of an event happening

	for (int i=0;i<N2;i++){
		if (communities[i]->nextAvailableTimeSlot == timeslot) communities[i]->isAvailable = true;
		if (communities[i]->isAvailable){
			int perturb = rand()%100;
			if (perturb == 0) randomlyPerturb(communities[i], timeslot, i);
		}
		if ((communities[i]->originFlag==2) && (communities[i]->birthTime >= timeslot - 100)){
			for (int j=0;j<10;j++) randomlyPerturb(communities[i], timeslot,i);
		}
	}

	if (z<=probEvent){
		generateEvent(timeslot);
	}
	cout << "Finished event here " << N2 << ", " << communities.size() << endl << flush;
	bool x = 0;
	for (int i = 0; i < N2; i++)  x &= communities[i]->isAvailable;
	cout << "After dummy check " << endl << flush;
	int numAC = 0;
	for (int i = 0; i < N2; i++){
		if (communities[i]->isAvailable) numAC++;
		cout << i << " " << flush;
	}

	cout << "Just a check " << communities[N2-1]->isAvailable << endl << flush;
	int numAvailableCommunities = 0;
	for (int i = 0; i < N2; i++){
		if (communities[i]->isAvailable) numAvailableCommunities++;
		cout << i << " " << flush;
	}
	cout << "Number of available communities at " << timeslot << " = " << numAvailableCommunities << endl << flush;
}
/********************************************************************/


void printGraph(){
	ofstream myfile;
	//myfile.open("ckbDynamicOut");
	myfile.open("ckbDynamicGraphByNode");
	for (int i=0;i<N1;i++){
		myfile << "Node: " << i << endl;
		for (int j=0;j<graph[i]->adj.size();j++){
			myfile << "(" << graph[i]->adj[j]->destId << ", "
				<< graph[i]->adj[j]->startTime << ", "
				<< ((graph[i]->adj[j]->endTime >= 0)? graph[i]->adj[j]->endTime : 2*T) << ", "
				<< graph[i]->adj[j]->communityId << ")" << endl;
			if (graph[i]->adj[j]->startTime <= -1 && graph[i]->adj[j]->communityId != -1) cout << "Yep, trouble " << graph[i]->adj[j]->startTime << ", " << graph[i]->adj[j]->communityId <<endl;
		}

	}
	myfile.close();
}

void printCommunity(){
	ofstream myfile;
//	myfile.open("ckbDynamicCommunitiesOut");
	myfile.open("ckbDynamicCoverByNode");
	/*int *numNodesTouched = new int[N2];
	for (int i=0;i<N2;i++) numNodesTouched[i] = (communodeInCommunitynities[i].nodeList).size();*/
	for (int i=0;i<N1;i++){
		if (graph[i]->communities.size() == 0) continue;
		myfile << i << ":";
		for (int j=0;j<graph[i]->communities.size();j++){
			int communityId = graph[i]->communities[j];
			community c = *communities[communityId];
			int pos = c.indexOfNode(i);
			if (pos==-1){
				cout << "Trouble with " << c.originFlag << endl << flush;
			}
			//numNodesTouched[communityId] -= 1;
			myfile << "(" << graph[i]->communities[j] << ", "
					<< c.nodeList[pos]->joinTime << ", "
					<< ((c.nodeList[pos]->leaveTime >= 0)? c.nodeList[pos]->leaveTime : 2*T)
				<< ")" << endl;
		}
		myfile << endl;
	}
	myfile.close();
	ofstream myfile1;
//	myfile1.open("ckbDynamicCommunitiesOut1");
	myfile1.open("ckbDynamicCoverByCommunity");
	/*int *numNodesTouched = new int[N2];
	for (int i=0;i<N2;i++) numNodesTouched[i] = (communities[i].nodeList).size();*/
	for (int i=0;i<N2;i++){
		myfile1 << i << ":";
		for (int j=0;j< (communities[i]->nodeList).size();j++){
			myfile1 << "(" << (communities[i]->nodeList[j])->nodeId
				<< ", " << (communities[i]->nodeList[j])->joinTime << ", "
				<< (((communities[i]->nodeList[j])->leaveTime >= 0)? (communities[i]->nodeList[j])->leaveTime : 2*T)
				<< ")" << endl;
		}
		myfile1 << endl;
	}
	myfile1.close();
}



void printGraphStream(){
	vector<update> stream;

	for (int i=0;i<N1;i++){
		if (graph[i]->startTime > 0)
			stream.push_back(update(1, i, -1, graph[i]->startTime));
		if (graph[i]->endTime > 0)
			stream.push_back(update(3, i, -1, graph[i]->endTime));

		for (int j=0;j<(graph[i]->adj).size();j++){
			if ((graph[i]->adj[j]->communityId == -1)||(graph[i]->adj[j]->communityId == -4))
				continue;
			if ((graph[i]->adj[j]->startTime > graph[i]->adj[j]->endTime) && (graph[i]->adj[j]->startTime != -1) && (graph[i]->adj[j]->endTime != -1)){
				cout << "OH MY GOD! THIS IS A PROBLEM! " << graph[i]->adj[j]->communityId << ", " << graph[i]->adj[j]->startTime << ", " << graph[i]->adj[j]->endTime << endl << flush;
				cout << "Edge: " << graph[i]->adj[j]->sourceId << ", " << graph[i]->adj[j]->destId << endl << flush;
				exit(1);
			}
			if (graph[i]->adj[j]->startTime == graph[i]->adj[j]->endTime) continue;
			if (graph[i]->adj[j]->startTime > 0)
				stream.push_back(update(2, i, graph[i]->adj[j]->destId, graph[i]->adj[j]->startTime));
			if (graph[i]->adj[j]->endTime > 0)
				stream.push_back(update(0, i, graph[i]->adj[j]->destId, graph[i]->adj[j]->endTime));
		}
	}
	sort(stream.begin(),stream.end(),compareUpdate);
	ofstream myfile;
	myfile.open("ckbDynamicStream");
	int first = 0, second = 1;
	while (second < stream.size()){
        update u1 = stream[first];
        update u2 = stream[second];
        if (u1.u == u2.u && u1.v == u2.v && u1.t == u2.t){
            first = second + 1;
            second = first + 1;
        }
        else{
            myfile << u1.updateType << "," << u1.u << "," << u1.v << "," << u1.t << endl;
            first = second;
            second = first + 1;
        }
	}
	if (first < stream.size())
        myfile << stream[first].updateType << "," << stream[first].u << "," << stream[first].v << "," << stream[first].t << endl;
	myfile.close();
}

void printInitialGraph(){
	ofstream myfile, myfileC;
//	myfile.open("ckbDynamicO");
//	myfileC.open("ckbDynamicCommunitiesO");
	myfile.open("ckbDynamicInitialGraphEdgeList");
	myfileC.open("ckbDynamicInitialCoverEdgeList");
	for (int i=0;i<N1;i++){
		for (int j=0; j < (graph[i]->adj).size(); j++){
			if (((graph[i]->adj[j])->communityId == -1)||((graph[i]->adj[j])->communityId == -4)) continue;
			if (((graph[i]->adj[j])->startTime <=0) && ((graph[i]->adj[j])->endTime >= -1 )){
				myfile << i << "\t" << graph[i]->adj[j]->destId << endl;
			}
		}

		for (int j=0;j<graph[i]->communities.size();j++){
			int communityId = graph[i]->communities[j];
			community c = *communities[communityId];
			int pos = c.indexOfNode(i);
			if (pos!=-1){
				if ((c.nodeList[pos]->joinTime >= 0) && (c.nodeList[pos]->leaveTime >= -1))
					myfileC << i << "\t" << communityId << endl;
			}
		}
	}
	myfile.close();
	myfileC.close();
}

void printAFOCSGraphAndStream(){
	ofstream myfile;
	myfile.open("AFOCS/tempGraph_0.txt");
	for (int i=0;i<N1;i++){
		for (int j=0; j < (graph[i]->adj).size(); j++){
			if (((graph[i]->adj[j])->communityId == -1)||((graph[i]->adj[j])->communityId == -4)) continue;
			if (((graph[i]->adj[j])->startTime <=0) && ((graph[i]->adj[j])->endTime >= -1 )){
				myfile << (i + 1) << " " << ((graph[i]->adj[j])->destId + 1) << endl;
			}
		}
	}
	myfile.close();
	for (int t = 0; t <= T; t++){
		std::ostringstream strs;
		strs << "AFOCS/tempGraph_";
		strs << (t+1);
		strs << ".txt";
		std::string outfilename = strs.str();
		//string outfilename =  to_string(alpha);
		const char* conv_outfilename = outfilename.c_str();
		myfile.open(conv_outfilename);
		for (int i = 0; i < N1; i++){
			for (int j = 0; j < (graph[i]->adj).size(); j++){
				if (((graph[i]->adj[j])->communityId == -1)||((graph[i]->adj[j])->communityId == -4)) continue;
				if ((graph[i]->adj[j])->startTime == t){
					myfile << (i + 1) << " " << ((graph[i]->adj[j])->destId + 1) << endl;
				}
				if ((graph[i]->adj[j])->endTime == t){
					myfile << "-" << (i + 1) << " " << ((graph[i]->adj[j])->destId + 1) << endl;
				}
			}
		}
		myfile.close();
	}
	myfile.open("AFOCS/tempGraphComm.txt");
	for (int i = 0; i < N2; i++){
		community c = *communities[i];
		bool flag = false;
		for (int j = 0; j < c.nodeList.size(); j++){
			if ((c.nodeList[j]->joinTime <= T) && (c.nodeList[j]->leaveTime == -1)){
				flag = true;
				myfile << (i + 1) << " ";
			}
			if (flag) myfile << endl;
		}
	}
	myfile.close();
}

int main(int argc, char *argv[]){
	if (argc < 9){
		cout << "Usage: ./ckbDynamicNodeSet numberOfNodes minCommunitySize maxCommunitySize minCommunityMembership ";
		cout << "maxCommunityMembership eventProbability intraCommunityEdgeProbability epsilon" << endl;
		exit(0);
	}
	N1 = atoi(argv[1]);
	mmin = atoi(argv[2]);
	mmax = atoi(argv[3]);
	xmin = atoi(argv[4]);
	xmax = atoi(argv[5]);
	probEvent = atof(argv[6]);
	alpha = atof(argv[7]);
	eps = atof(argv[8]);


	srand(time(NULL));
	cout << "Initializing Sampler " << N1 << ", " << xmin << ", " << xmax << ", " << beta1 << endl << flush;
	sampler = new BucketSampling(0, xmin, xmax, beta1, rand());

	cout << "Initialized Sampler " << endl << flush;
	struct timespec start, finish, mid;
	double elapsedS = 0, elapsedD = 0;
	clock_gettime(CLOCK_MONOTONIC, &start);
	cout << "Generating static structure ... " << endl << flush;
	generateStaticStructure();
	clock_gettime(CLOCK_MONOTONIC, &mid);
	cout << "Generated static structure. " << endl << flush;
	for (int t = 0; t< T;t++){
		generateTimeslot(t);
		/*bool isok = sanityCheck2();
		if (!isok){
            cout << "Sanity check 2 is false at " << t << endl;
            exit(0);
		}*/
		/*Number of nodes with 1 and 2 memberships, and total memberships*/
		int n1 = 0, n2 = 0, total = 0;
		for (int i = 0; i < N1; i++){
			if (nodeMemberships[i] == 1) n1 += 1;
			if (nodeMemberships[i] == 2) n2 += 1;
			total += nodeMemberships[i];
		}
		assert(total == currMem);
		double x = initialMem/currMem;
		probB = 0.5*x/(1+x);	//linear
		//probB = (0.5 + x*x)/(1 + x*x); //quadratic
		probD = 0.5 - probB;
		cout << "currMem: " << currMem << " initialMem: " << initialMem << " probB: " << probB << " probD: " << probD << endl;
		cout << "orphan : " << t << " " << numOrphanNodes << " " << nBirths << " " << nDeaths
			<< " " << n1 << " " << n2 << " " << total << endl << flush;
	}
	clock_gettime(CLOCK_MONOTONIC, &finish);
	elapsedS = (mid.tv_sec - start.tv_sec) + (mid.tv_nsec - start.tv_nsec)/pow(10,9);
	elapsedD = (finish.tv_sec - start.tv_sec) + (finish.tv_nsec - start.tv_nsec)/pow(10,9);
	printGraph();
	printCommunity();

	printInitialGraph();
	printGraphStream();
	printAFOCSGraphAndStream();
	averageBirthTime /= nBirths;
	averageDeathTime /= nDeaths;
	averageSplitTime /= nSplits;
	averageMergeTime /= nMerges;
	averageNodeAddTime /= nAdd;
	averageNodeDeleteTime /= nDelete;
	cout << N1 << "\t" << probEvent << "\t" << alpha << "\t" << eps << "\t" << elapsedS << "\t" << elapsedD << endl;
}

