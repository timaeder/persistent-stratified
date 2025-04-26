#pragma once
#include <iostream>
#include <list>
#include <vector>
#include <set>
#include <cmath>
#include <algorithm>
#include "EuclDist.hh"

/**
 * @brief Represents a neighborhood graph for a dataset.
 * 
 * @tparam T Numeric type for the dataset and distances (e.g., float, double, long double).
 */
template <typename T>
class NeighborhoodGraph {
private:
    std::vector<std::vector<T>> dataset;  // Input dataset (points in space)
    T epsilon;                            // Neighborhood radius

    // Adjacency graphs
    std::vector<std::vector<int>> graph1; // Graph with edges for j >= i
    std::vector<std::vector<int>> graph2; // Graph with all edges

public:
    /**
     * @brief Constructor for the NeighborhoodGraph class.
     * @param dataset Input dataset as a vector of vectors (points in space).
     * @param epsilon Neighborhood radius.
     */
    NeighborhoodGraph(const std::vector<std::vector<T>>& dataset, const T& epsilon)
        : dataset(dataset), epsilon(epsilon) {}

    /**
     * @brief Computes the standard neighborhood graph.
     */
    void computeStandard() {
        size_t n = dataset.size();
        graph1.resize(n);
        graph2.resize(n);

        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                T dist = vectorDistance(dataset[i].begin(), dataset[i].end(), dataset[j].begin());
                if (dist <= epsilon) {
                    if (j >= i) {
                        graph1[i].push_back(j);
                    }
                    graph2[i].push_back(j);
                }
            }
        }
    }

    /**
     * @brief Computes the CVR (Constrained Voronoi Region) neighborhood graph.
     */
    void computeCVR() {
        size_t k = dataset.size();
        std::vector<std::set<int>> tempGraph1(k + 1);
        std::vector<std::set<int>> tempGraph2(k + 1);

        for (size_t i = 0; i < k; i++) {
            tempGraph2[i].insert(i);
            tempGraph1[i].insert(i);

            for (size_t j = i + 1; j < k; j++) {
                T dist = 0.5 * vectorDistance(dataset[i].begin(), dataset[i].end(), dataset[j].begin());
                if (dist < 0.5 * epsilon) {
                    tempGraph1[i].insert(j);
                    tempGraph2[i].insert(j);
                    tempGraph2[j].insert(i);
                }
            }
            tempGraph1[i].insert(k);
            tempGraph2[i].insert(k);
            tempGraph2[k].insert(i);
        }

        tempGraph2[k].insert(k);
        tempGraph1[k].insert(k);

        graph1.resize(tempGraph1.size());
        graph2.resize(tempGraph2.size());

        for (size_t i = 0; i < tempGraph1.size(); i++) {
            graph1[i] = std::vector<int>(tempGraph1[i].begin(), tempGraph1[i].end());
            graph2[i] = std::vector<int>(tempGraph2[i].begin(), tempGraph2[i].end());
        }
    }

    /**
     * @brief Computes the boundary-projected neighborhood graph.
     */
    void computeBoundary() {
        size_t k = dataset.size();
        std::vector<T> m(dataset[0].size());
        std::vector<std::set<int>> tempGraph1(k);
        std::vector<std::set<int>> tempGraph2(k);

        for (size_t i = 0; i < k; i++) {
            if ((0.5 * epsilon - vectorDistance(dataset[i].begin(), dataset[i].end(), dataset[0].begin())) < 0.5 * epsilon) {
                tempGraph2[i].insert(i);
                tempGraph1[i].insert(i);

                for (size_t j = i + 1; j < k; j++) {
                    // Compute midpoint (m)
                    typename std::vector<T>::const_iterator first = dataset[i].begin();
                    typename std::vector<T>::const_iterator first2 = dataset[j].begin();

                    for (size_t it = 0; it < m.size(); it++) {
                        m[it] = ((*first++) + (*first2++)) / 2;
                    }

                    T dist = (vectorDistance(m.begin(), m.end(), dataset[0].begin()) < 0.5 * epsilon)
                                 ? BoundaryProjectedDistance(dataset[i].begin(), dataset[i].end(), dataset[j].begin(), dataset[0].begin(), dataset[0].end(), 0.5 * epsilon)
                                 : 0.5 * vectorDistance(dataset[i].begin(), dataset[i].end(), dataset[j].begin());

                    if (dist < 0.5 * epsilon) {
                        tempGraph1[i].insert(j);
                        tempGraph2[i].insert(j);
                        tempGraph2[j].insert(i);
                    }
                }
            }
        }

        graph1.resize(tempGraph1.size());
        graph2.resize(tempGraph2.size());

        for (size_t i = 0; i < tempGraph1.size(); i++) {
            graph1[i] = std::vector<int>(tempGraph1[i].begin(), tempGraph1[i].end());
            graph2[i] = std::vector<int>(tempGraph2[i].begin(), tempGraph2[i].end());
        }
    }

    /**
     * @brief Gets the adjacency graphs.
     * @return A pair of adjacency graphs (graph1, graph2).
     */
    std::pair<std::vector<std::vector<int>>, std::vector<std::vector<int>>> getGraphs() const {
        return {graph1, graph2};
    }
};
