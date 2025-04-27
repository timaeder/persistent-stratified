#pragma once
#include <vector>
#include <list>
#include <set>
#include <algorithm>
#include <utility>
#include "my_hasher.hh"
#include "DistanceCalculator.hh"
#include "newNeighborhood.hh"
#include "ksubset.hh"
#include "unordered_map.hpp"

class Complex {
public:
    std::vector<std::vector<std::pair<long double, std::vector<int>>>> simplices;
    std::list<std::vector<std::vector<int>>> neighborhood;
    ska::unordered_map<std::vector<int>, int, VectorHasher> IDs;

    // Compute the diameter of a simplex (standard)
    long double computeSimplexDiameter(const std::vector<std::vector<long double>>& points, const std::vector<int>& simplex) {
        int d = simplex.size();
        long double maxDiameter = 0.0;

        if (d > 1) {
            for (int i = 0; i < d; i++) {
                for (int j = i + 1; j < d; j++) {
                    long double dist = DistanceCalculator::vectorDistance(
                        points[simplex[i]].begin(), points[simplex[i]].end(),
                        points[simplex[j]].begin()
                    );
                    maxDiameter = std::max(maxDiameter, dist);
                }
            }
        }
        return maxDiameter;
    }

    // Compute the diameter of a simplex (relative)
    long double computeRelSimplexDiameter(const std::vector<std::vector<long double>>& points, const std::vector<int>& simplex, const long double epsilon) {
        int d = simplex.size();
        long double maxDiameter = 0.0;

        if (d > 1) {
            for (int i = 0; i < d; i++) {
                for (int j = i + 1; j < d; j++) {
                    long double dist = DistanceCalculator::vectorDistance(
                        points[simplex[i]].begin(), points[simplex[i]].end(),
                        points[simplex[j]].begin()
                    );
                    maxDiameter = std::max(maxDiameter, dist);
                }
            }
        }

        // Ensure the diameter is less than or equal to 0.5 * epsilon
        return maxDiameter < 0.5 * epsilon ? maxDiameter : epsilon;
    }

    // Compute the diameter of a simplex (boundary)
    long double computeBndSimplexDiameter(const std::vector<std::vector<long double>>& points, const std::vector<int>& simplex, const long double epsilon) {
        int d = simplex.size();
        long double maxDiameter = 0.0;
        std::vector<long double> midpoint(points[0].size());

        if (d > 1) {
            for (int i = 0; i < d; i++) {
                for (int j = i + 1; j < d; j++) {
                    // Compute the midpoint
                    auto first = points[simplex[i]].begin();
                    auto last = points[simplex[i]].end();
                    auto first2 = points[simplex[j]].begin();

                    for (size_t it = 0; it < midpoint.size(); it++) {
                        midpoint[it] = ((*first++) + (*first2++)) / 2;
                    }

                    // Check if the midpoint is within the boundary
                    long double dist;
                    if (DistanceCalculator::vectorDistance(midpoint.begin(), midpoint.end(), points[0].begin()) < 0.5 * epsilon) {
                        dist = DistanceCalculator::boundaryProjectedDistance(
                            points[simplex[i]].begin(), points[simplex[i]].end(),
                            points[simplex[j]].begin(),
                            points[0].begin(), points[0].end(),
                            0.5 * epsilon
                        );
                    } else {
                        dist = 0.5 * DistanceCalculator::vectorDistance(points[simplex[i]].begin(), points[simplex[i]].end(), points[simplex[j]].begin());
                    }

                    maxDiameter = std::max(maxDiameter, dist);
                }
            }
        } else {
            long double temp = 0.5 * epsilon - DistanceCalculator::vectorDistance(points[simplex[0]].begin(), points[simplex[0]].end(), points[0].begin());
            if (temp > 0) {
                maxDiameter = temp;
            }
        }

        return maxDiameter;
    }

    // Add cofaces for IncrementalVR
    void addCofaces(
        const std::vector<std::vector<long double>>& points,
        const std::vector<std::vector<int>>& graph,
        int k,
        const std::vector<int>& simplex,
        const std::vector<int>& neighbors
    ) {
        int d = simplex.size();
        simplices[d - 1].push_back({computeSimplexDiameter(points, simplex), simplex});

        if (d > k + 1) return;

        for (const int& neighbor : neighbors) {
            if (std::find(simplex.begin(), simplex.end(), neighbor) == simplex.end()) {
                std::vector<int> newSimplex = simplex;
                newSimplex.push_back(neighbor);
                std::sort(newSimplex.begin(), newSimplex.end());

                std::vector<int> intersection;
                std::set_intersection(
                    neighbors.begin(), neighbors.end(),
                    graph[neighbor].begin(), graph[neighbor].end(),
                    std::back_inserter(intersection)
                );

                addCofaces(points, graph, k, newSimplex, intersection);
            }
        }
    }

    // Add cofaces for RelIncrementalVR
    void addRelCofaces(
        const std::vector<std::vector<long double>>& points,
        const std::vector<std::vector<int>>& graph,
        int k,
        const std::vector<int>& simplex,
        const std::vector<int>& neighbors,
        const long double epsilon
    ) {
        int d = simplex.size();
        long double diameter = computeRelSimplexDiameter(points, simplex, epsilon);

        if (diameter < 0.5 * epsilon) {
            simplices[d - 1].push_back({diameter, simplex});
        } else {
            return;
        }

        if (d < k + 2) {
            for (const int& neighbor : neighbors) {
                if (std::find(simplex.begin(), simplex.end(), neighbor) == simplex.end()) {
                    std::vector<int> newSimplex = simplex;
                    newSimplex.push_back(neighbor);
                    std::sort(newSimplex.begin(), newSimplex.end());

                    std::vector<int> intersection;
                    std::set_intersection(
                        neighbors.begin(), neighbors.end(),
                        graph[neighbor].begin(), graph[neighbor].end(),
                        std::back_inserter(intersection)
                    );

                    addRelCofaces(points, graph, k, newSimplex, intersection, epsilon);
                }
            }
        }
    }

    // Add cofaces for BndIncrementalVR
    void addBndCofaces(
        const std::vector<std::vector<long double>>& points,
        const std::vector<std::vector<int>>& graph,
        int k,
        const std::vector<int>& simplex,
        const std::vector<int>& neighbors,
        const long double epsilon
    ) {
        int d = simplex.size();
        simplices[d - 1].push_back({computeBndSimplexDiameter(points, simplex, epsilon), simplex});

        if (d < k + 2) {
            for (const int& neighbor : neighbors) {
                if (std::find(simplex.begin(), simplex.end(), neighbor) == simplex.end()) {
                    std::vector<int> newSimplex = simplex;
                    newSimplex.push_back(neighbor);
                    std::sort(newSimplex.begin(), newSimplex.end());

                    std::vector<int> intersection;
                    std::set_intersection(
                        neighbors.begin(), neighbors.end(),
                        graph[neighbor].begin(), graph[neighbor].end(),
                        std::back_inserter(intersection)
                    );

                    addBndCofaces(points, graph, k, newSimplex, intersection, epsilon);
                }
            }
        }
    }

    // Constructor
    Complex(int dimension) : simplices(dimension) {}

    // Set the neighborhood
    void setNeighborhood(const std::list<std::vector<std::vector<int>>>& nbh) {
        neighborhood = nbh;
    }

    // Compute IncrementalVR
    void computeIncrementalVR(
        const std::vector<std::vector<long double>>& points,
        const long double& epsilon,
        const int& k
    ) {
        Neighborhood neighborhoodCalculator(points, epsilon);
        setNeighborhood(neighborhoodCalculator.computeNeighborhood());

        int n = neighborhood.front().size();
        std::vector<int> simplex(1);

        for (int i = 0; i < n; i++) {
            simplex[0] = i;
            addCofaces(points, neighborhood.front(), k, simplex, neighborhood.front()[i]);
        }

        int count = 0;
        for (int i = k + 1; i >= 0; --i) {
            std::stable_sort(
                simplices[i].begin(), simplices[i].end(),
                [](const auto& left, const auto& right) { return left.first < right.first; }
            );

            for (auto it = simplices[i].rbegin(); it != simplices[i].rend(); ++it) {
                IDs[it->second] = count++;
            }
        }
    }

    // Compute RelIncrementalVR
    void computeRelIncrementalVR(
        const std::vector<std::vector<long double>>& points,
        const long double& epsilon,
        const int& k
    ) {
        Neighborhood neighborhoodCalculator(points, epsilon);
        setNeighborhood(neighborhoodCalculator.computeNeighborhood());

        int n = neighborhood.front().size();
        std::vector<int> simplex(1);

        for (int i = 0; i < n; i++) {
            if (!neighborhood.front()[i].empty()) {
                simplex[0] = neighborhood.front()[i][0];
                addRelCofaces(points, neighborhood.front(), k, simplex, neighborhood.front()[simplex[0]], epsilon);
            }
        }

        int count = 0;
        for (int i = k + 1; i >= 0; --i) {
            std::stable_sort(
                simplices[i].begin(), simplices[i].end(),
                [](const auto& left, const auto& right) { return left.first < right.first; }
            );

            for (auto it = simplices[i].rbegin(); it != simplices[i].rend(); ++it) {
                IDs[it->second] = count++;
            }
        }
    }

    // Compute BndIncrementalVR
    void computeBndIncrementalVR(
        const std::vector<std::vector<long double>>& points,
        const long double& epsilon,
        const int& k
    ) {
        Neighborhood neighborhoodCalculator(points, epsilon);
        setNeighborhood(neighborhoodCalculator.computeNeighborhood());

        int n = neighborhood.front().size();
        std::vector<int> simplex(1);

        for (int i = 0; i < n; i++) {
            if (!neighborhood.front()[i].empty()) {
                simplex[0] = neighborhood.front()[i][0];
                addBndCofaces(points, neighborhood.front(), k, simplex, neighborhood.front()[simplex[0]], epsilon);
            }
        }

        int count = 0;
        for (int i = k + 1; i >= 0; --i) {
            std::stable_sort(
                simplices[i].begin(), simplices[i].end(),
                [](const auto& left, const auto& right) { return left.first < right.first; }
            );

            for (auto it = simplices[i].rbegin(); it != simplices[i].rend(); ++it) {
                IDs[it->second] = count++;
            }
        }
    }
};
