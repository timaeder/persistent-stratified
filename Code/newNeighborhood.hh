#pragma once
#include <vector>
#include <list>
#include <set>
#include "DistanceCalculator.hh"

class Neighborhood {
private:
    std::vector<std::vector<long double>> points;
    long double epsilon;

public:
    // Constructor
    Neighborhood(const std::vector<std::vector<long double>>& points, long double epsilon)
        : points(points), epsilon(epsilon) {}

    // Compute the standard neighborhood
    std::list<std::vector<std::vector<int>>> computeNeighborhood() {
        int n = points.size();
        std::vector<std::vector<int>> graph1(n);
        std::vector<std::vector<int>> graph2(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                long double dist = DistanceCalculator::vectorDistance(points[i].begin(), points[i].end(), points[j].begin());
                if (dist <= epsilon) {
                    if (j >= i) {
                        graph1[i].push_back(j);
                    }
                    graph2[i].push_back(j);
                }
            }
        }

        std::list<std::vector<std::vector<int>>> graphs;
        graphs.push_back(graph1);
        graphs.push_back(graph2);
        return graphs;
    }

    // Compute the CVR neighborhood
    std::list<std::vector<std::vector<int>>> computeCVRNeighborhood() {
        int k = points.size();
        std::vector<std::set<int>> graph1(k + 1);
        std::vector<std::set<int>> graph2(k + 1);

        for (int i = 0; i < k; i++) {
            graph2[i].insert(i);
            graph1[i].insert(i);

            for (int j = i + 1; j < k; j++) {
                long double dist = 0.5 * DistanceCalculator::vectorDistance(points[i].begin(), points[i].end(), points[j].begin());
                if (dist < 0.5 * epsilon) {
                    graph1[i].insert(j);
                    graph2[i].insert(j);
                    graph2[j].insert(i);
                }
            }
            graph1[i].insert(k);
            graph2[i].insert(k);
            graph2[k].insert(i);
        }

        graph2[k].insert(k);
        graph1[k].insert(k);

        std::vector<std::vector<int>> graph3;
        std::vector<std::vector<int>> graph4;

        for (const auto& g : graph1) {
            graph3.push_back(std::vector<int>(g.begin(), g.end()));
        }

        for (const auto& g : graph2) {
            graph4.push_back(std::vector<int>(g.begin(), g.end()));
        }

        std::list<std::vector<std::vector<int>>> graphs;
        graphs.push_back(graph3);
        graphs.push_back(graph4);
        return graphs;
    }

    // Compute the boundary neighborhood
    std::list<std::vector<std::vector<int>>> computeBndNeighborhood() {
        int k = points.size();
        std::vector<long double> m(points[0].size());
        std::vector<std::set<int>> graph1(k);
        std::vector<std::set<int>> graph2(k);

        for (int i = 0; i < k; i++) {
            if ((0.5 * epsilon - DistanceCalculator::vectorDistance(points[i].begin(), points[i].end(), points[0].begin())) < 0.5 * epsilon) {
                graph2[i].insert(i);
                graph1[i].insert(i);

                for (int j = i + 1; j < k; j++) {
                    // Compute midpoint (m)
                    auto first = points[i].begin();
                    auto last = points[i].end();
                    auto first2 = points[j].begin();

                    for (size_t it = 0; it < m.size(); it++) {
                        m[it] = ((*first++) + (*first2++)) / 2;
                    }

                    long double dist;
                    if (DistanceCalculator::vectorDistance(m.begin(), m.end(), points[0].begin()) < 0.5 * epsilon) {
                        dist = DistanceCalculator::boundaryProjectedDistance(
                            points[i].begin(), points[i].end(),
                            points[j].begin(),
                            points[0].begin(), points[0].end(),
                            0.5 * epsilon
                        );
                    } else {
                        dist = 0.5 * DistanceCalculator::vectorDistance(points[i].begin(), points[i].end(), points[j].begin());
                    }

                    if (dist < 0.5 * epsilon) {
                        graph1[i].insert(j);
                        graph2[i].insert(j);
                        graph2[j].insert(i);
                    }
                }
            }
        }

        // Convert graph1 and graph2 from sets to vectors
        std::vector<std::vector<int>> graph3;
        std::vector<std::vector<int>> graph4;

        for (const auto& g : graph1) {
            graph3.push_back(std::vector<int>(g.begin(), g.end()));
        }

        for (const auto& g : graph2) {
            graph4.push_back(std::vector<int>(g.begin(), g.end()));
        }

        std::list<std::vector<std::vector<int>>> graphs;
        graphs.push_back(graph3);
        graphs.push_back(graph4);
        return graphs;
    }
};
