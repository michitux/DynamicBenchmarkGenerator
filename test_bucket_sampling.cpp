#include "BucketSampling.h"
#include <map>

void printDistribution(const std::vector<int>& values) {
	std::map<int, int> hist_map;

	for (int deg : values) {
		++hist_map[deg];
	}

	for (const auto& it : hist_map) {
		std::cout << "deg " << it.first << ": " << it.second << std::endl;
	}
}

int main() {
	int n = 10000;
	BucketSampling sampler(0, 1, n/10, -2.5, 42);

	PowerlawDegreeSequence memberships_sequence(1, n/10, -2.5);
	memberships_sequence.run();

	std::vector<int> desiredMemberships(memberships_sequence.getDegreeSequence(n));

	for (int i = 0; i < n; ++i) {
		sampler.addNode(desiredMemberships[i]);
	}

	std::vector<int> actualMemberships(sampler.getNumberOfNodes());
	int actualSumOfMemberships = 0;

	PowerlawDegreeSequence community_size_sequence(std::min(20, n/100), n/10, -2.5);
	community_size_sequence.run();

	while (actualSumOfMemberships < sampler.getSumOfDesiredMemberships()) {
		int wanted_size = community_size_sequence.getDegree();
		std::vector<int> community = sampler.birthCommunityNodes(wanted_size);

		actualSumOfMemberships += community.size();

		assert(community.size() == wanted_size);

		for (int u : community) {
			sampler.assignCommunity(u);
			++actualMemberships[u];
		}
	}

	std::cout << "For " << n << " wanted nodes we have " << sampler.getNumberOfNodes() << " actual nodes (due to oversampling)" << std::endl;

	printDistribution(desiredMemberships);
	printDistribution(actualMemberships);

	// try removing all nodes again
	for (int i = 0; i < sampler.getNumberOfNodes(); ++i) {
		if (sampler.getDesiredMemberships(i) > 0) {
			sampler.removeNode(i);
		}
	}

	return 0;
};
