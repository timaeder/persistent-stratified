#pragma once
#include <iostream>
#include <vector>
#include <list>
#include <unordered_map>
#include <set>
#include <cmath>
#include <algorithm>
#include "unordered_map.hpp"
#include "Neighborhood.hh"
#include "my_hasher.hh"
#include "ksubset.hh"

/**
 * @brief Represents a simplex in a simplicial complex.
 */
class Simplex {
public:
    std::vector<int> vertices;  // Indices of the vertices forming the simplex
    long double diameter;       // Diameter of the simplex

    /**
     * @brief Constructor for a simplex.
     * @param vertices Indices of the vertices forming the simplex.
     * @param diameter Diameter of the simplex.
     */
    Simplex(const std::vector<int>& vertices, long double diameter)
        : vertices(vertices), diameter(diameter) {}

    /**
     * @brief Computes the diameter of a simplex.
     * @param S Input dataset as a vector of vectors (points in space).
     * @return The diameter of the simplex.
     */
    template <typename T>
    static T computeDiameter(const std::vector<std::vector<T>>& S, const std::vector<int>& vertices) {
        size_t d = vertices.size();
        T maxDiameter = T(0);

        if (d > 1) {
            for (size_t i = 0; i < d; i++) {
                for (size_t j = i + 1; j < d; j++) {
                    T dist = vectorDistance(S[vertices[i]].begin(), S[vertices[i]].end(), S[vertices[j]].begin());
                    if (dist > maxDiameter) {
                        maxDiameter = dist;
                    }
                }
            }
        }
        return maxDiameter;
    }
};

/**
 * @brief Represents a simplicial complex.
 */
class SimplicialComplex {
public:
    std::vector<std::vector<Simplex>> simplices;  // Simplices organized by dimension
    ska::unordered_map<std::vector<int>, int, VectorHasher> IDs;  // Map of simplex IDs

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
     * @brief Constructs an incremental Vietoris-Rips complex.
     * @param S Input dataset as a vector of vectors (points in space).
     * @param epsilon Neighborhood radius.
     * @param k Maximum dimension of the simplicial complex.
     * @param neighborhoodGraph NeighborhoodGraph object.
     * @return The constructed Vietoris-Rips complex.
     */
    template <typename T>
    static SimplicialComplex constructIncrementalVR(
        const std::vector<std::vector<T>>& S,
        const T& epsilon,
        const int& k,
        NeighborhoodGraph<T>& neighborhoodGraph
    ) {
        // Compute the standard neighborhood graph
        neighborhoodGraph.computeStandard();
        auto [graph1, graph2] = neighborhoodGraph.getGraphs();

        SimplicialComplex VR(k + 1);

        size_t n = graph1.size();
        std::vector<int> u(1);

        for (size_t i = 0; i < n; i++) {
            u[0] = i;
            addCofaces(S, graph1, k, u, graph1[i], VR);
        }

        VR.finalize();
        return VR;
    }

private:
    /**
     * @brief Adds cofaces to the simplicial complex.
     * @param S Input dataset as a vector of vectors (points in space).
     * @param G Neighborhood graph (adjacency list).
     * @param k Maximum dimension of the simplicial complex.
     * @param t Current simplex.
     * @param N Neighboring points.
     * @param complex Simplicial complex to update.
     */
    template <typename T>
    static void addCofaces(
        const std::vector<std::vector<T>>& S,
        const std::vector<std::vector<int>>& G,
        const int& k,
        const std::vector<int>& t,
        const std::vector<int>& N,
        SimplicialComplex& complex
    ) {
        int d = t.size();
        complex.addSimplex(Simplex(t, Simplex::computeDiameter(S, t)), d - 1);

        if (d > k + 1) {
            return;
        }

        for (const int& f : N) {
            if (std::find(t.begin(), t.end(), f) == t.end()) {
                std::vector<int> M;
                std::vector<int> s = t;
                s.push_back(f);
                std::set_intersection(
                    N.begin(), N.end(),
                    G[f].begin(), G[f].end(),
                    std::back_inserter(M)
                );
                addCofaces(S, G, k, s, M, complex);
            }
        }
    }
};
