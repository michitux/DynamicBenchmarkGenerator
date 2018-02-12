#pragma once

#include <cstddef>
#include <cassert>
#include <vector>
#include <array>
#include <iostream>
#include <random>
#include <algorithm>

#include "PowerlawDegreeSequence.h"

class BucketSampling {
private:
	static constexpr int oversample_fraction = 5;

	struct node_data {
		// positions in the full_slots vector where this node occurs
		std::vector<int> full_slot_positions;
		// fractional part of the number of slots times oversample_fraction
		// this is negative if the node is in too many communities
		int current_fractional_slot;
		// position in the fractional slot vector of current_fractional_slot > 0
		int fractional_slot_position;
		// desired degree of the node. 0 if the node has been removed.
		int degree;
	};

	struct slot {
		int node; /* The node for which the slot is */
		int pos_at_node; /* The position in the node's slot position vector that links to the slot */
	};

	PowerlawDegreeSequence pl_generator;
	// distribution of the actually generated degrees
	std::vector<int> degree_distribution;
	std::vector<std::array<int, oversample_fraction-1>> nodes_with_fractional_slots_for_degree;
	std::vector<node_data> nodes;
	std::vector<bool> node_sampled_in_current_call;
	std::vector<slot> full_slots;
	std::array<std::vector<int>, oversample_fraction> fractional_slots;
	std::mt19937_64 random_engine;

	int sumOfDesiredMemberships;

	void verifyInvariants(bool verify_oversampling = false) {
#ifndef NDEBUG
		int num_full_slots = 0;
		std::array<int, oversample_fraction> num_fractional_slots;
		num_fractional_slots.fill(0);

		int num_existing_nodes = 0;

		for (int u = 0; u < nodes.size(); ++u) {
			node_data& node = nodes[u];

			if (node.degree == 0) {
				assert(node.full_slot_positions.empty());
				assert(node.current_fractional_slot <= 0);
			} else {
				++num_existing_nodes;

				if (node.current_fractional_slot < 0) {
					assert(node.full_slot_positions.empty());
				}

				for (int i = 0; i < node.full_slot_positions.size(); ++i) {
					const int slot_pos = node.full_slot_positions[i];
					++num_full_slots;

					assert(full_slots[slot_pos].node == u);
					assert(full_slots[slot_pos].pos_at_node == i);
				}

				if (node.current_fractional_slot > 0) {
					++num_fractional_slots[node.current_fractional_slot];
					assert(fractional_slots[node.current_fractional_slot][node.fractional_slot_position] == u);

					auto& nodes_with_frac = nodes_with_fractional_slots_for_degree[node.degree];
					assert(std::find(nodes_with_frac.begin(), nodes_with_frac.end(), u) != nodes_with_frac.end());
				}
			}
		}

		assert(num_full_slots == full_slots.size());
		for (int i = 0; i < oversample_fraction; ++i) {
			assert(num_fractional_slots[i] == fractional_slots[i].size());
		}

		int membership_sum = 0;
		for (int i = 0; i <= pl_generator.getMaximumDegree(); ++i) {
			if (i < pl_generator.getMinimumDegree()) {
				assert(degree_distribution[i] == 0);
				for (int v : nodes_with_fractional_slots_for_degree[i]) {
					assert(v == -1);
				}
			} else if (verify_oversampling) {
				int nodes_found = 0;

				for (int v : nodes_with_fractional_slots_for_degree[i]) {
					if (v != -1) {
						assert(nodes[v].degree == i);
						++nodes_found;
					}
				}

				assert(nodes_found == (degree_distribution[i] % (oversample_fraction + 1)));
			}

			membership_sum += degree_distribution[i];
		}

		assert(membership_sum == num_existing_nodes);
#endif
	}
public:
	BucketSampling(int n, int minSlots, int maxSlots, double exponent, uint64_t seed) : pl_generator(minSlots, maxSlots, exponent), degree_distribution(maxSlots + 1, 0), nodes_with_fractional_slots_for_degree(maxSlots + 1), random_engine(seed), sumOfDesiredMemberships(0) {
		pl_generator.run();

		for (int d = 0; d <= maxSlots; ++d) {
			nodes_with_fractional_slots_for_degree[d].fill(-1);
		}

		for (int u = 0; u < n; ++u) {
			addNode();
		}
	}

	std::vector<int> birthCommunityNodes(int communitySize) {
		std::vector<int> result;

		std::array<double, oversample_fraction> slot_boundaries;

		slot_boundaries[0] = full_slots.size();
		for (int i = 1; i < oversample_fraction; ++i) {
			// The weight of each slot in this array
			double slot_weight = static_cast<double>(i) / oversample_fraction;
			double total_slot_weight = static_cast<double>(fractional_slots[i].size()) * slot_weight;
			slot_boundaries[i] = slot_boundaries[i-1] + total_slot_weight;
		}

		const double max_slot = slot_boundaries.back();
		std::uniform_real_distribution<double> distr(0.0, max_slot);

		while (result.size() < communitySize) {
			const double sampled_slot = distr(random_engine);
			int u = 0;

			if (sampled_slot < slot_boundaries[0]) {
				u = full_slots[static_cast<int>(sampled_slot)].node;
			} else {
				for (int i = 1; i < oversample_fraction; ++i) {
					if (sampled_slot < slot_boundaries[i]) {
						double stretch_factor = (static_cast<double>(oversample_fraction) / i);
						double fractional_pos_in_bin = sampled_slot - slot_boundaries[i - 1];
						u = fractional_slots[i][static_cast<int>(fractional_pos_in_bin * stretch_factor)];
						break;
					}
				}
			}

			if (!node_sampled_in_current_call[u]) {
				result.push_back(u);
				node_sampled_in_current_call[u] = true;
			}
		}

		for (int u : result) {
			node_sampled_in_current_call[u] = false;
		}

		return result;
	}

	void assignCommunity(int nodeId) {
		node_data& node = nodes[nodeId];

		if (!node.full_slot_positions.empty()) {
			// relocate slot to the end of full_slots
			int slot_pos = node.full_slot_positions.back();
			assert(full_slots[slot_pos].node == nodeId);

			std::swap(full_slots[slot_pos], full_slots.back());

			// repair the reference to the other slot position
			slot& other_slot = full_slots[slot_pos];
			nodes[other_slot.node].full_slot_positions[other_slot.pos_at_node] = slot_pos;
			node.full_slot_positions.pop_back();
			full_slots.pop_back();
		} else {
			removeFraction(nodeId, oversample_fraction);
		}

		verifyInvariants();
	}

	void leaveCommunity(int nodeId) {
		node_data& node = nodes[nodeId];

		/* The node was removed internally but was still in a community. Do not give it back any slots. */
		if (node.degree == 0) return;

		if (node.current_fractional_slot < 0) {
			assert(node.full_slot_positions.empty());

			addFraction(nodeId, oversample_fraction);
		} else {
			int slot_pos = node.full_slot_positions.size();
			node.full_slot_positions.push_back(full_slots.size());
			full_slots.push_back({nodeId, slot_pos});
		}

		verifyInvariants();
	}

private:
	void removeFraction(int nodeId, int fraction) {
		if (fraction == 0) return;

		node_data& node = nodes[nodeId];

		if (node.current_fractional_slot > 0) {
			std::vector<int>& frac_slots = fractional_slots[node.current_fractional_slot];
			// relocate slot to the end of the fractional slots
			int slot_pos = node.fractional_slot_position;
			std::swap(frac_slots[slot_pos], frac_slots.back());

			int other_node = frac_slots[slot_pos];
			// repair the reference to the other slot position
			nodes[other_node].fractional_slot_position = slot_pos;
			// remove the slot
			frac_slots.pop_back();
		}

		node.current_fractional_slot -= fraction;
	}

	void addFraction(int nodeId, int fraction) {
		if (fraction == 0) return;

		node_data& node = nodes[nodeId];

		node.current_fractional_slot += fraction;

		if (node.current_fractional_slot > 0) {
			std::vector<int>& frac_slots = fractional_slots[node.current_fractional_slot];
			node.fractional_slot_position = frac_slots.size();
			frac_slots.push_back(nodeId);
		}
	}

	int addNode(int degree, bool oversample) {
		const int u = nodes.size();
		nodes.emplace_back();
		node_sampled_in_current_call.push_back(false);

		node_data& node = nodes[u];

		node.degree = degree;
		++degree_distribution[node.degree];

		int u_fractional_slot = node.degree % oversample_fraction;
		int u_additional_full_slots = node.degree / oversample_fraction;
		int num_full_slots = node.degree;

		if (oversample) {
			auto& nodes_with_fractional_slots = nodes_with_fractional_slots_for_degree[node.degree];
			if (degree_distribution[node.degree] % (oversample_fraction + 1) == oversample_fraction) {
				// remove oversampled slots
				for (int& v : nodes_with_fractional_slots) {
					assert(v != -1);

					for (int i = 0; i < u_additional_full_slots; ++i) {
						assignCommunity(v);
					}

					removeFraction(v, u_fractional_slot);
					v = -1;
				}

				// create a new node with exactly degree slots
				addNode(degree, false);

				// no fractional slot for this node
				u_fractional_slot = 0;
			} else {
				num_full_slots += u_additional_full_slots;
				addFraction(u, u_fractional_slot);

				bool found = false;
				for (int& v : nodes_with_fractional_slots) {
					if (v == -1) {
						v = u;
						found = true;
						break;
					}
				}

				assert(found);
			}
		} else {
			assert(degree_distribution[node.degree] % (oversample_fraction + 1) == 0);
			u_fractional_slot = 0;
		}

		node.full_slot_positions.reserve(num_full_slots);

		for (int i = 0; i < num_full_slots; ++i) {
			leaveCommunity(u);
		}


		assert(node.full_slot_positions.size() == num_full_slots);
		assert(node.current_fractional_slot == u_fractional_slot);

		return u;
	}

public:

	int addNode() {
		return addNode(pl_generator.getDegree());
	}

	int addNode(int degree) {
		int u = addNode(degree, true);
		sumOfDesiredMemberships += degree;
		verifyInvariants(true);
		return u;
	}

private:
	void removeNode(int nodeId, bool withFraction) {
		node_data& node = nodes[nodeId];

		// make sure there are no more slots for this node
		while (node.current_fractional_slot >= 0) {
			assignCommunity(nodeId);
		}

		const int degree = node.degree;

		--degree_distribution[degree];
		node.degree = 0;

		auto& nodes_with_fractional_slots = nodes_with_fractional_slots_for_degree[degree];
		bool hadFraction = false;
		int nodesWithFraction = 0;
		for (int& v : nodes_with_fractional_slots) {
			if (v == nodeId) {
				v = -1;
				hadFraction = true;
			} else if (v != -1) {
				++nodesWithFraction;
			}
		}

		assert(withFraction || !hadFraction);

		// if node had a fraction, we are
		// Otherwise, we need to get the fractional part from somewhere else.
		if (!hadFraction && withFraction) {
			int u_fractional_slot = degree % oversample_fraction;
			int u_additional_full_slots = degree / oversample_fraction;

			if (nodesWithFraction > 0) {
				// if there is another node with a fractional part, it is easy, too:
				// just remove the additional slots from this node.
				std::uniform_int_distribution<int> random_selector(0, nodesWithFraction-1);
				int pos = random_selector(random_engine);

				int i = 0;
				for (int& v : nodes_with_fractional_slots) {
					if (v != -1) {
						if (i == pos) {
							// remove additional fractional part from this node
							for (int j = 0; j < u_additional_full_slots; ++j) {
								assignCommunity(j);
							}

							removeFraction(v, u_fractional_slot);

							v = -1;
							break;
						}
						++i;
					}
				}
			} else {
				// now the difficult case: we need to
				// a) remove a complete other node and
				// b) find oversample_fraction nodes of the same degree that then get a fraction!

				// FIXME store nodes of degree for every degree!
				std::vector<int> nodes_of_degree;

				for (int u = 0; u < nodes.size(); ++u) {
					if (nodes[u].degree == degree) {
						nodes_of_degree.push_back(u);
					}
				}

				auto pop_random_node = [&]() -> int {
					std::uniform_int_distribution<int> sample(0, nodes_of_degree.size() - 1);
					int pos = sample(random_engine);
					std::swap(nodes_of_degree[pos], nodes_of_degree.back());
					int x = nodes_of_degree.back();
					nodes_of_degree.pop_back();
					return x;
				};

				assert(nodes_of_degree.size() >= oversample_fraction);

				verifyInvariants();

				removeNode(pop_random_node(), false);

				verifyInvariants();

				// add a fraction to oversample_fraction nodes of nodes_of_degree
				for (int& u : nodes_with_fractional_slots) {
					assert(u == -1);

					u = pop_random_node();

					assert(nodes[u].degree == degree);

					for (int j = 0; j < u_additional_full_slots; ++j) {
						leaveCommunity(u);
					}

					addFraction(u, u_fractional_slot);
				}
			}
		}
	}
public:
	void removeNode(int nodeId) {
		if (nodes[nodeId].degree == 0) return;

		sumOfDesiredMemberships -= nodes[nodeId].degree;
		removeNode(nodeId, true);
		verifyInvariants(true);
	}

	int getNumberOfNodes() const {
		return nodes.size();
	}

	int getSumOfDesiredMemberships() const {
		return sumOfDesiredMemberships;
	}

	int getDesiredMemberships(int nodeId) const {
		return nodes[nodeId].degree;
	}
};
