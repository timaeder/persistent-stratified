#pragma once
#include <iostream>
#include <fstream>
#include <list>
#include <unordered_map>
#include <set>
#include <vector>
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
     * @param vertices Vertices of the simplex.
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

    /**
     * @brief Computes the boundary-projected diameter of a simplex.
     * @param S Input dataset as a vector of vectors (points in space).
     * @param vertices Vertices of the simplex.
     * @param epsilon Neighborhood radius.
     * @return The boundary-projected diameter of the simplex.
     */
    static long double computeBoundarySimplexDiameter(
        const std::vector<std::vector<long double>>& S,
        const std::vector<int>& vertices,
        const long double& epsilon
    ) {
        int d = vertices.size();
        long double maxDiameter = 0.0;
        std::vector<long double> m(S[0].size());

        if (d > 1) {
            for (int i = 0; i < d; i++) {
                for (int j = i + 1; j < d; j++) {
                    auto first = S[vertices[i]].begin();
                    auto first2 = S[vertices[j]].begin();
                    for (size_t it = 0; it < m.size(); it++) {
                        m[it] = ((*first++) + (*first2++)) / 2;
                    }

                    long double dist;
                    if (vectorDistance(m.begin(), m.end(), S[0].begin()) < 0.5 * epsilon) {
                        dist = BoundaryProjectedDistance(S[vertices[i]].begin(), S[vertices[i]].end(), S[vertices[j]].begin(), S[0].begin(), S[0].end(), 0.5 * epsilon);
                    } else {
                        dist = 0.5 * vectorDistance(S[vertices[i]].begin(), S[vertices[i]].end(), S[vertices[j]].begin());
                    }
                    maxDiameter = std::max(maxDiameter, dist);
                }
            }
        } else {
            long double dist = 0.5 * epsilon - vectorDistance(S[vertices[0]].begin(), S[vertices[0]].end(), S[0].begin());
            if (dist > 0) {
                maxDiameter = dist;
            }
        }
        return maxDiameter;
    }

    /**
     * @brief Computes the relative diameter of a simplex.
     * @param S Input dataset as a vector of vectors (points in space).
     * @param vertices Vertices of the simplex.
     * @param epsilon Neighborhood radius.
     * @return The relative diameter of the simplex.
     */
    static long double computeRelativeSimplexDiameter(
        const std::vector<std::vector<long double>>& S,
        const std::vector<int>& vertices,
        const long double& epsilon
    ) {
        int d = vertices.size();
        long double maxDiameter = 0.0;
        std::vector<std::vector<int>> T;
        std::vector<int> t;
        std::vector<int> s2;
        std::vector<long double> m(S[0].size());

        if (d > 2) {
            if (vertices[d - 1] >= S.size()) {
                s2 = vertices;
                s2.erase(s2.end() - 1);
                subset(s2, d - 1, 2, 0, t, T);

                for (const auto& sigma : T) {
                    auto first = S[sigma[0]].begin();
                    auto first2 = S[sigma[1]].begin();
                    for (int it = 0; it < m.size(); it++) {
                        m[it] = ((*first++) + (*first2++)) / 2;
                    }

                    long double dist;
                    if (vectorDistance(m.begin(), m.end(), S[0].begin()) < 0.5 * epsilon) {
                        dist = BoundaryProjectedDistance(S[sigma[0]].begin(), S[sigma[0]].end(), S[sigma[1]].begin(), S[0].begin(), S[0].end(), 0.5 * epsilon);
                    } else {
                        dist = 0.5 * vectorDistance(S[sigma[0]].begin(), S[sigma[0]].end(), S[sigma[1]].begin());
                    }
                    maxDiameter = std::max(maxDiameter, dist);
                }
            } else {
                subset(vertices, d, 2, 0, t, T);
                for (const auto& sigma : T) {
                    long double dist = 0.5 * vectorDistance(S[sigma[0]].begin(), S[sigma[0]].end(), S[sigma[1]].begin());
                    maxDiameter = std::max(maxDiameter, dist);
                }
            }
        } else if (d > 1) {
            if (vertices[1] < S.size()) {
                maxDiameter = 0.5 * vectorDistance(S[vertices[0]].begin(), S[vertices[0]].end(), S[vertices[1]].begin());
            } else {
                long double dist = (0.5 * epsilon - vectorDistance(S[0].begin(), S[0].end(), S[vertices[0]].begin()));
                maxDiameter = std::max(static_cast<long double>(0.0), dist);
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
    std::list<std::vector<std::vector<int>>> neighborhood; // Neighborhood structure

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
     * @brief Adds cofaces to the complex.
     * @param S Input dataset.
     * @param G Neighborhood graph.
     * @param k Maximum dimension.
     * @param t Current simplex.
     * @param N Neighboring vertices.
     */
    void addCofaces(
        const std::vector<std::vector<long double>>& S,
        const std::vector<std::vector<int>>& G,
        const int& k,
        const std::vector<int>& t,
        const std::vector<int>& N
    ) {
        int d = t.size();
        simplices[d - 1].emplace_back(Simplex::computeDiameter(S, t), t);

        if (d > k + 1) {
            return;
        }

        for (const int& f : N) {
            if (std::find(t.begin(), t.end(), f) == t.end()) {
                std::vector<int> M;
                std::vector<int> s = t;
                s.push_back(f);
                std::set_intersection(N.begin(), N.end(), G[f].begin(), G[f].end(), std::back_inserter(M));
                addCofaces(S, G, k, s, M);
            }
        }
    }

    /**
     * @brief Adds relative cofaces to the complex.
     * @param S Input dataset.
     * @param G Neighborhood graph.
     * @param k Maximum dimension.
     * @param t Current simplex.
     * @param N Neighboring vertices.
     * @param epsilon Neighborhood radius.
     */
    void addRelativeCofaces(
        const std::vector<std::vector<long double>>& S,
        const std::vector<std::vector<int>>& G,
        const int& k,
        const std::vector<int>& t,
        const std::vector<int>& N,
        const long double epsilon
    ) {
        int d = t.size();
        long double diam = Simplex::computeRelativeSimplexDiameter(S, t, epsilon);

        if (diam < 0.5 * epsilon) {
            simplices[d - 1].emplace_back(diam, t);
        } else {
            return;
        }

        if (d < k + 2) {
            for (const int& f : N) {
                if (std::find(t.begin(), t.end(), f) == t.end()) {
                    std::vector<int> M;
                    std::vector<int> s = t;
                    s.push_back(f);
                    std::sort(s.begin(), s.end());

                    if (s.back() < S.size()) {
                        std::set_intersection(N.begin(), N.end(), G[f].begin(), G[f].end(), std::back_inserter(M));
                        addRelativeCofaces(S, G, k, s, M, epsilon);
                    } else {
                        addRelativeCofaces(S, G, k, s, M, epsilon);
                    }
                }
            }
        }
    }

    /**
     * @brief Adds boundary cofaces to the complex.
     * @param S Input dataset.
     * @param G Neighborhood graph.
     * @param k Maximum dimension.
     * @param t Current simplex.
     * @param N Neighboring vertices.
     * @param epsilon Neighborhood radius.
     */
    void addBoundaryCofaces(
        const std::vector<std::vector<long double>>& S,
        const std::vector<std::vector<int>>& G,
        const int& k,
        const std::vector<int>& t,
        const std::vector<int>& N,
        const long double epsilon
    ) {
        int d = t.size();
        simplices[d - 1].emplace_back(Simplex::computeBoundarySimplexDiameter(S, t, epsilon), t);

        if (d < k + 2) {
            for (const int& f : N) {
                if (std::find(t.begin(), t.end(), f) == t.end()) {
                    std::vector<int> M;
                    std::vector<int> s = t;
                    s.push_back(f);
                    std::sort(s.begin(), s.end());
                    std::set_intersection(N.begin(), N.end(), G[f].begin(), G[f].end(), std::back_inserter(M));
                    addBoundaryCofaces(S, G, k, s, M, epsilon);
                }
            }
        }
    }

    /**
     * @brief Constructs an incremental Vietoris-Rips complex.
     * @param S Input dataset.
     * @param epsilon Neighborhood radius.
     * @param k Maximum dimension of the simplicial complex.
     * @param nbh Neighborhood graph.
     * @return The constructed Vietoris-Rips complex.
     */
    SimplicialComplex constructIncrementalVR(
        const std::vector<std::vector<long double>>& S,
        const long double epsilon,
        const int k,
        const std::list<std::vector<std::vector<int>>>& nbh
    ) {
        SimplicialComplex VR(k + 2);
        VR.neighborhood = nbh;

        int n = VR.neighborhood.front().size();
        std::vector<int> u(1);

        for (int i = 0; i < n; i++) {
            u[0] = i;
            addCofaces(S, VR.neighborhood.front(), k, u, VR.neighborhood.front()[i]);
        }

        VR.finalize();
        return VR;
    }

    /**
     * @brief Constructs a relative incremental Vietoris-Rips complex.
     * @param S Input dataset.
     * @param k Maximum dimension of the simplicial complex.
     * @param nbh Neighborhood graph.
     * @param epsilon Neighborhood radius.
     * @return The constructed relative Vietoris-Rips complex.
     */
    SimplicialComplex constructRelativeIncrementalVR(
        const std::vector<std::vector<long double>>& S,
        const int k,
        const std::list<std::vector<std::vector<int>>>& nbh,
        const long double epsilon
    ) {
        SimplicialComplex VR(k + 2);
        VR.neighborhood = nbh;

        int n = VR.neighborhood.front().size();
        std::vector<int> u(1);

        for (int i = 0; i < n; i++) {
            if (!VR.neighborhood.front()[i].empty()) {
                u[0] = VR.neighborhood.front()[i][0];
                addRelativeCofaces(S, VR.neighborhood.front(), k, u, VR.neighborhood.front()[u[0]], epsilon);
            }
        }

        VR.finalize();
        return VR;
    }

    /**
     * @brief Constructs a boundary incremental Vietoris-Rips complex.
     * @param S Input dataset.
     * @param k Maximum dimension of the simplicial complex.
     * @param nbh Neighborhood graph.
     * @param epsilon Neighborhood radius.
     * @return The constructed boundary Vietoris-Rips complex.
     */
    SimplicialComplex constructBoundaryIncrementalVR(
        const std::vector<std::vector<long double>>& S,
        const int k,
        const std::list<std::vector<std::vector<int>>>& nbh,
        const long double epsilon
    ) {
        SimplicialComplex VR(k + 2);
        VR.neighborhood = nbh;

        int n = VR.neighborhood.front().size();
        std::vector<int> u(1);

        for (int i = 0; i < n; i++) {
            if (!VR.neighborhood.front()[i].empty()) {
                u[0] = VR.neighborhood.front()[i][0];
                addBoundaryCofaces(S, VR.neighborhood.front(), k, u, VR.neighborhood.front()[u[0]], epsilon);
            }
        }

        VR.finalize();
        return VR;
    }

    /**
     * @brief Computes coboundary of a simplex.
     * @param s Simplex vertices.
     * @param N Neighborhood graph.
     * @param n Maximum dimension.
     * @return Coboundary indices.
     */
    std::vector<int> computeCoboundary(
        const std::vector<int>& s,
        const std::vector<std::vector<int>>& N,
        const int& n
    ) {
        int k = s.size() - 1;
        std::vector<int> cobound;
        std::vector<int> last_intersection = N[s[0]];
        std::vector<int> curr_intersection;

        if (last_intersection.size() > k + 1) {
            bool stop = false;
            for (int i = 1; i < s.size(); i++) {
                if (N[s[i]].size() > 0 && !stop) {
                    std::set_intersection(last_intersection.begin(), last_intersection.end(),
                        N[s[i]].begin(), N[s[i]].end(),
                        std::back_inserter(curr_intersection));
                    std::swap(last_intersection, curr_intersection);
                    curr_intersection.clear();
                } else {
                    stop = true;
                    last_intersection.clear();
                }
            }
            if (k < n) {
                if (last_intersection.size()) {
                    for (const int& v : last_intersection) {
                        if (std::find(s.begin(), s.end(), v) == s.end()) {
                            std::vector<int> t = s;
                            t.push_back(v);
                            std::sort(t.begin(), t.end());
                            cobound.push_back(IDs.find(t)->second);
                        }
                    }
                }
            }
            if (cobound.size()) {
                std::sort(cobound.begin(), cobound.end());
            }
        }
        return cobound;
    }

    /**
     * @brief Computes relative coboundary of a simplex.
     * @param S Input dataset.
     * @param s Simplex vertices.
     * @param N Neighborhood graph.
     * @param n Maximum dimension.
     * @param epsilon Neighborhood radius.
     * @return Relative coboundary indices.
     */
    std::vector<int> computeRelativeCoboundary(
        const std::vector<std::vector<long double>>& S,
        const std::vector<int>& s,
        const std::vector<std::vector<int>>& N,
        const int& n,
        const long double epsilon
    ) {
        int k = s.size() - 1;
        std::vector<int> cobound;
        std::vector<int> last_intersection = N[s[0]];
        std::vector<int> curr_intersection;
        int conept = S.size();

        if (last_intersection.size() > k + 1) {
            bool stop = false;
            for (int i = 1; i < s.size(); i++) {
                if (N[s[i]].size() > 0 && !stop) {
                    std::set_intersection(last_intersection.begin(), last_intersection.end(),
                        N[s[i]].begin(), N[s[i]].end(),
                        std::back_inserter(curr_intersection));
                    std::swap(last_intersection, curr_intersection);
                    curr_intersection.clear();
                } else {
                    stop = true;
                    last_intersection.clear();
                }
            }
            if (k < n) {
                if (last_intersection.size()) {
                    for (const int& v : last_intersection) {
                        if (std::find(s.begin(), s.end(), v) == s.end()) {
                            std::vector<int> t = s;
                            t.push_back(v);
                            std::sort(t.begin(), t.end());
                            if (t.back() != conept) {
                                cobound.push_back(IDs.find(t)->second);
                            } else {
                                if (Simplex::computeDiameter(S, t) < 0.5 * epsilon) {
                                    cobound.push_back(IDs.find(t)->second);
                                }
                            }
                        }
                    }
                }
            }
            if (cobound.size()) {
                std::sort(cobound.begin(), cobound.end());
            }
        }
        return cobound;
    }

    /**
     * @brief Adds boundary to a simplex.
     * @param cbID1 First coboundary ID vector.
     * @param cbID2 Second coboundary ID vector.
     */
    void addBoundary(std::vector<int>& cbID1, std::vector<int>& cbID2) {
        std::vector<int> symm_diff;
        std::set_symmetric_difference(
            cbID1.begin(), cbID1.end(),
            cbID2.begin(), cbID2.end(),
            std::back_inserter(symm_diff)
        );
        std::swap(cbID1, symm_diff);
    }

    /**
     * @brief Computes copairings for a complex.
     * @param S Input dataset.
     * @param epsilon Neighborhood radius.
     * @param k Maximum dimension of the simplicial complex.
     * @param nbh Neighborhood graph.
     * @return Vector of copairings.
     */
    std::vector<std::vector<std::pair<long double, long double>>> computeCopairings(
        const std::vector<std::vector<long double>>& S,
        const long double& epsilon,
        const int& k,
        const std::list<std::vector<std::vector<int>>>& nbh
    ) {
        int n = simplices.size();
        std::vector<std::vector<std::pair<long double, long double>>> copairs(n);

        int full_size = IDs.size();
        int curr_size = 0;

        std::vector<std::vector<int>> redcols(full_size, std::vector<int>{-1});

        for (int d = 0; d < n - 1; d++) {
            std::vector<std::set<int>> redrows(full_size);
            curr_size += simplices[d].size();

            for (auto simp_it1 = simplices[d].rbegin(); simp_it1 != simplices[d].rend(); ++simp_it1) {
                const std::vector<int>& simp1 = simp_it1->vertices;
                int id1 = full_size - curr_size + (simp_it1 - simplices[d].rbegin());

                if (redcols[id1].size()) {
                    auto cbID1 = computeCoboundary(simp1, nbh.back(), n - 1);
                    redcols[id1] = cbID1;

                    if (cbID1.size()) {
                        int pivot1 = cbID1.back();
                        bool stop = false;

                        while (redrows[pivot1].size() > 0 && !stop) {
                            auto cbID2 = redcols[*redrows[pivot1].begin()];
                            addBoundary(cbID1, cbID2);

                            if (cbID1.size()) {
                                pivot1 = cbID1.back();
                            } else {
                                stop = true;
                            }
                        }
                    }

                    redcols[id1] = cbID1;

                    if (cbID1.size()) {
                        int pivot1 = cbID1.back();  // Declare and initialize pivot1 here
                        redrows[pivot1].insert(id1);
                        if (redcols[*cbID1.rbegin()].size()) {
                            if (std::abs(simplices[d + 1][full_size - curr_size - *cbID1.rbegin() - 1].diameter - simp_it1->diameter) > 0.0001) {
                                copairs[d].emplace_back(simp_it1->diameter, simplices[d + 1][full_size - curr_size - *cbID1.rbegin() - 1].diameter);
                            }
                            redcols[*cbID1.rbegin()].clear();
                        }
                    } else {
                        copairs[d].emplace_back(simp_it1->diameter, epsilon + 0.1);
                    }
                }
            }
        }

        return copairs;
    }

    /**
     * @brief Computes local copairings for a complex.
     * @param S Input dataset.
     * @param k Maximum dimension of the simplicial complex.
     * @param z Central point.
     * @param rad Radius.
     * @return Vector of local copairings.
     */
    std::vector<std::vector<std::pair<long double, long double>>> computeLocalCopairings(
        const std::vector<std::vector<long double>>& S,
        const int& k,
        const std::vector<long double>& z,
        const long double rad
    ) {
        std::vector<std::pair<long double, std::vector<long double>>> FiltPts;
        long double epsilon = rad;

        for (const auto& v : S) {
            long double r = vectorDistance(z.begin(), z.end(), v.begin());
            if (r <= rad) {
                FiltPts.emplace_back(r, v);
            }
        }

        std::stable_sort(FiltPts.begin(), FiltPts.end());

        std::vector<std::vector<long double>> LPC;
        for (const auto& v : FiltPts) {
            LPC.push_back(v.second);
        }

        // Create an instance of the NeighborhoodGraph class with the dataset and epsilon
        NeighborhoodGraph<long double> neighborhoodGraph(LPC, rad);

        // Compute the CVR neighborhood graph
        neighborhoodGraph.computeCVR();

        // Retrieve the computed graphs
        std::list<std::vector<std::vector<int>>> nbh(neighborhoodGraph.getGraphs().first.begin(), neighborhoodGraph.getGraphs().first.end()); // or .second based on your requirement

        // Use the computed neighborhood in your VR construction
        SimplicialComplex VR = constructRelativeIncrementalVR(LPC, k, nbh, rad);

        int n = k + 1;
        for (int i = 0; i <= k + 1; i++) {
            if (!VR.simplices[i].size()) {
                n = i - 1;
                break;
            }
        }

        int full_size = VR.IDs.size();
        int curr_size = 0;

        std::vector<std::vector<int>> redcols(full_size, std::vector<int>{-1});
        std::vector<std::set<int>> redrows(full_size);

        std::vector<std::vector<std::pair<long double, long double>>> LocCopairs(k + 1);

        for (int d = 0; d < n; d++) {
            curr_size += VR.simplices[d].size();

            for (auto simp_it1 = VR.simplices[d].rbegin(); simp_it1 != VR.simplices[d].rend(); ++simp_it1) {
                const std::vector<int>& simp1 = simp_it1->vertices;
                int id1 = full_size - curr_size + (simp_it1 - VR.simplices[d].rbegin());

                if (redcols[id1].size()) {
                    auto cbID1 = computeRelativeCoboundary(LPC, simp1, VR.neighborhood.back(), n, rad);
                    redcols[id1] = cbID1;

                    if (cbID1.size()) {
                        int pivot1 = cbID1.back();
                        bool stop = false;

                        while (redrows[pivot1].size() > 0 && !stop) {
                            auto cbID2 = redcols[*redrows[pivot1].begin()];
                            addBoundary(cbID1, cbID2);

                            if (cbID1.size()) {
                                pivot1 = cbID1.back();
                            } else {
                                stop = true;
                            }
                        }
                    }

                    redcols[id1] = cbID1;

                    if (cbID1.size()) {
                        int pivot1 = cbID1.back();  // Declare and initialize pivot1 here
                        redrows[pivot1].insert(id1);
                        if (redcols[*cbID1.rbegin()].size()) {
                            long double pers = std::fabs(VR.simplices[d + 1][full_size - curr_size - *cbID1.rbegin() - 1].diameter - simp_it1->diameter);
                            if (pers > 0.001) {
                                LocCopairs[d].emplace_back(simp_it1->diameter, pers);
                            }
                            redcols[*cbID1.rbegin()].clear();
                        }
                    } else {
                        LocCopairs[d].emplace_back(simp_it1->diameter, (0.5 * rad) - simp_it1->diameter + 0.001);
                    }
                }
            }
        }

        return LocCopairs;
    }

    /**
     * @brief Computes boundary-projected local copairings for a complex.
     * @param S Input dataset.
     * @param k Maximum dimension of the simplicial complex.
     * @param z Central point.
     * @param rad Radius.
     * @return Vector of boundary-projected local copairings.
     */
    std::vector<std::vector<std::pair<long double, long double>>> computeBoundaryLocalCopairings(
        const std::vector<std::vector<long double>>& S,
        const int& k,
        const std::vector<long double>& z,
        const long double rad
    ) {
        std::vector<std::pair<long double, std::vector<long double>>> FiltPts;
        long double epsilon = rad;

        for (const auto& v : S) {
            long double r = vectorDistance(z.begin(), z.end(), v.begin());
            if (r <= rad) {
                FiltPts.emplace_back(r, v);
            }
        }

        std::stable_sort(FiltPts.begin(), FiltPts.end());

        std::vector<std::vector<long double>> LPC;
        for (const auto& v : FiltPts) {
            LPC.push_back(v.second);
        }

        // Create an instance of the NeighborhoodGraph class for the dataset and radius
        NeighborhoodGraph<long double> neighborhoodGraph(LPC, rad);

        // Compute the boundary-projected neighborhood graph
        neighborhoodGraph.computeBoundary();

        // Convert the resulting vector to a list
        std::list<std::vector<std::vector<int>>> nbh(neighborhoodGraph.getGraphs().first.begin(), neighborhoodGraph.getGraphs().first.end());

        // Construct the boundary incremental VR complex using the converted neighborhood graph
        SimplicialComplex VR = constructBoundaryIncrementalVR(LPC, k, nbh, rad);

        int n = k + 1;
        for (int i = 0; i <= k + 1; i++) {
            if (!VR.simplices[i].size()) {
                n = i - 1;
                break;
            }
        }

        int full_size = VR.IDs.size();
        int curr_size = 0;
        std::vector<std::vector<int>> redcols(full_size, std::vector<int>{-1});
        std::vector<std::set<int>> redrows(full_size);
        std::vector<std::vector<std::pair<long double, long double>>> LocCopairs(k + 1);

        for (int d = 0; d < n; d++) {
            curr_size += VR.simplices[d].size();
            for (auto simp_it1 = VR.simplices[d].rbegin(); simp_it1 != VR.simplices[d].rend(); ++simp_it1) {
                const std::vector<int>& simp1 = simp_it1->vertices;
                int id1 = full_size - curr_size + (simp_it1 - VR.simplices[d].rbegin());

                if (redcols[id1].size()) {
                    auto cbID1 = computeCoboundary(simp1, VR.neighborhood.back(), n);
                    redcols[id1] = cbID1;

                    if (cbID1.size()) {
                        int pivot1 = cbID1.back();
                        bool stop = false;

                        while (redrows[pivot1].size() > 0 && !stop) {
                            auto cbID2 = redcols[*redrows[pivot1].begin()];
                            addBoundary(cbID1, cbID2);

                            if (cbID1.size()) {
                                pivot1 = cbID1.back();
                            } else {
                                stop = true;
                            }
                        }
                    }

                    redcols[id1] = cbID1;

                    if (cbID1.size()) {
                        int pivot1 = cbID1.back();  // Declare and initialize pivot1 here
                        redrows[pivot1].insert(id1);
                    
                        if (redcols[*cbID1.rbegin()].size()) {
                            long double pers = std::abs(VR.simplices[d + 1][full_size - curr_size - *cbID1.rbegin() - 1].diameter - simp_it1->diameter);
                            if (pers > 0.001) {
                                LocCopairs[d].emplace_back(simp_it1->diameter, pers);
                            }
                            redcols[*cbID1.rbegin()].clear();
                        }
                    } else {
                        LocCopairs[d].emplace_back(simp_it1->diameter, (0.5 * rad) - simp_it1->diameter + 0.001);
                    }
                }
            }
        }

        // "reduce homology"
        if (!LocCopairs[0].empty()) {
            LocCopairs[0].erase(LocCopairs[0].end() - 1);
        }

        return LocCopairs;
    }
};
