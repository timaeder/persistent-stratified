#pragma once
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <cmath>
#include <algorithm>
#include "RipsComplex.hh"

/**
 * @brief Represents a simplicial complex with methods for persistent cohomology.
 */
class SimplicialComplex {
public:
    std::vector<std::vector<Simplex>> simplices;  // Simplices organized by dimension
    ska::unordered_map<std::vector<int>, int, VectorHasher> IDs;  // Map of simplex IDs
    std::list<std::vector<std::vector<int>>> neighborhood;  // Neighborhood graphs

    /**
     * @brief Constructor for a simplicial complex.
     * @param maxDimension Maximum dimension of the simplicial complex.
     */
    SimplicialComplex(int maxDimension)
        : simplices(maxDimension + 1) {}

    /**
     * @brief Adds a simplex to the complex.
     * @param simplex The simplex to add.
     * @param dimension The dimension of the simplex.
     */
    void addSimplex(const Simplex& simplex, int dimension) {
        simplices[dimension].push_back(simplex);
    }

    /**
     * @brief Sorts simplices by diameter and assigns IDs.
     */
    void finalize() {
        int count = 0;
        for (int i = simplices.size() - 1; i >= 0; --i) {
            std::stable_sort(simplices[i].begin(), simplices[i].end(),
                [](const Simplex& left, const Simplex& right) {
                    return left.diameter < right.diameter;
                });

            for (const auto& simplex : simplices[i]) {
                IDs[simplex.vertices] = count++;
            }
        }
    }

    /**
     * @brief Computes the coboundary of a simplex.
     * @param simplexVertices Vertices of the simplex.
     * @param dimension Dimension of the simplex.
     * @return A vector of coboundary simplex IDs.
     */
    std::vector<int> computeCoboundary(const std::vector<int>& simplexVertices, int dimension) const {
        int k = simplexVertices.size() - 1;
        std::vector<int> coboundary;
        std::vector<int> lastIntersection = neighborhood.back()[simplexVertices[0]];
        std::vector<int> currentIntersection;

        if (lastIntersection.size() > k + 1) {
            for (size_t i = 1; i < simplexVertices.size(); i++) {
                std::set_intersection(
                    lastIntersection.begin(), lastIntersection.end(),
                    neighborhood.back()[simplexVertices[i]].begin(),
                    neighborhood.back()[simplexVertices[i]].end(),
                    std::back_inserter(currentIntersection)
                );
                std::swap(lastIntersection, currentIntersection);
                currentIntersection.clear();
            }

            if (k < dimension) {
                for (const int& vertex : lastIntersection) {
                    if (std::find(simplexVertices.begin(), simplexVertices.end(), vertex) == simplexVertices.end()) {
                        std::vector<int> newSimplex = simplexVertices;
                        newSimplex.push_back(vertex);
                        std::sort(newSimplex.begin(), newSimplex.end());
                        coboundary.push_back(IDs.at(newSimplex));
                    }
                }
            }
        }

        std::sort(coboundary.begin(), coboundary.end());
        return coboundary;
    }

    /**
     * @brief Computes persistent cohomology for the simplicial complex.
     * @param maxDimension Maximum dimension for cohomology computation.
     * @return A vector of persistence pairs for each dimension.
     */
    std::vector<std::vector<std::pair<long double, long double>>> computePersistentCohomology(int maxDimension) {
        int fullSize = IDs.size();
        std::vector<std::vector<int>> reductionColumns(fullSize, std::vector<int>{-1});
        std::vector<std::set<int>> reductionRows(fullSize);

        std::vector<std::vector<std::pair<long double, long double>>> persistencePairs(maxDimension + 1);

        for (int d = 0; d < maxDimension; d++) {
            int currentSize = 0;
            for (int i = 0; i <= d; i++) {
                currentSize += simplices[i].size();
            }

            for (auto simplexIt = simplices[d].rbegin(); simplexIt != simplices[d].rend(); ++simplexIt) {
                const std::vector<int>& simplexVertices = simplexIt->vertices;
                int simplexID = fullSize - currentSize + (simplexIt - simplices[d].rbegin());

                if (!reductionColumns[simplexID].empty()) {
                    std::vector<int> coboundary = computeCoboundary(simplexVertices, maxDimension);

                    if (!coboundary.empty()) {
                        int pivot = coboundary.back();
                        bool stop = false;

                        while (!reductionRows[pivot].empty() && !stop) {
                            std::vector<int> otherColumn = reductionColumns[*reductionRows[pivot].begin()];
                            std::vector<int> symmetricDifference;

                            std::set_symmetric_difference(
                                coboundary.begin(), coboundary.end(),
                                otherColumn.begin(), otherColumn.end(),
                                std::back_inserter(symmetricDifference)
                            );

                            coboundary = symmetricDifference;

                            if (!coboundary.empty()) {
                                pivot = coboundary.back();
                            } else {
                                stop = true;
                            }
                        }

                        reductionColumns[simplexID] = coboundary;

                        if (!coboundary.empty()) {
                            reductionRows[pivot].insert(simplexID);
                            if (!reductionColumns[coboundary.back()].empty()) {
                                long double birth = simplexIt->diameter;
                                long double death = simplices[d + 1][fullSize - currentSize - coboundary.back() - 1].diameter;

                                if (std::abs(death - birth) > 0.001) {
                                    persistencePairs[d].emplace_back(birth, death);
                                }

                                reductionColumns[coboundary.back()].clear();
                            }
                        } else {
                            persistencePairs[d].emplace_back(simplexIt->diameter, std::numeric_limits<long double>::infinity());
                        }
                    }
                }
            }
        }

        return persistencePairs;
    }
};
