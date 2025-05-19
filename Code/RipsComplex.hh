#pragma once
#include <iostream>
#include <fstream>
#include <list>
#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include "unordered_map.hpp"
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include "Neighborhood.hh"
#include "my_hasher.hh"
#include "ksubset.hh"


template <typename T>
struct Cplx {
    Cplx() : simplices(), neighborhood(), IDs() {}
    std::vector<std::vector<std::pair<T, std::vector<int>>>> simplices;
    std::list<std::vector<std::vector<int>>> neighborhood;
    ska::unordered_map<std::vector<int>, int, VectorHasher> IDs;

    Cplx(int n) {
        simplices = std::vector<std::vector<std::pair<T, std::vector<int>>>>(n);
    }
};

template <typename T>
T SimplexDiameter(const std::vector<std::vector<T>>& S, const std::vector<int>& s) {
    int d = s.size();
    T maxsofar = T{0};
    T temp = T{0};

    if (d > 1) {
        for (int i = 0; i < d; i++) {
            for (int j = i + 1; j < d; j++) {
                temp = vectorDistance(S[s[i]].begin(), S[s[i]].end(), S[s[j]].begin());
                if (temp > maxsofar) {
                    maxsofar = temp;
                }
            }
        }
    }
    return maxsofar;
}


template <typename T>
T BndSimplexDiameter(
    const std::vector<std::vector<T>>& S,
    const std::vector<int>& s,
    const T& epsilon
) {
    int d = s.size();
    T maxsofar = T{0};
    T temp = T{0};
    std::vector<T> m(S[0].size());

    if (d > 1) {
        for (int i = 0; i < d; i++) {
            for (int j = i + 1; j < d; j++) {
                // Compute midpoint (m)
                auto first = S[s[i]].begin();
                auto first2 = S[s[j]].begin();

                for (size_t it = 0; it < m.size(); it++) {
                    m[it] = ((*first++) + (*first2++)) / 2;
                }

                // Use the updated G function
                if (vectorDistance(m.begin(), m.end(), S[0].begin()) < T{0.5} * epsilon) {
                    temp = BoundaryProjectedDistance(S[s[i]].begin(), S[s[i]].end(), S[s[j]].begin(), S[0].begin(), S[0].end(), T{0.5} * epsilon);
                } else {
                    temp = T{0.5} * vectorDistance(S[s[i]].begin(), S[s[i]].end(), S[s[j]].begin());
                }

                if (temp > maxsofar) {
                    maxsofar = temp;
                }
            }
        }
    } else {
        temp = T{0.5} * epsilon - vectorDistance(S[s[0]].begin(), S[s[0]].end(), S[0].begin());
        if (temp > T{0}) {
            maxsofar = temp;
        }
    }
    return maxsofar;
}


template <typename T>
T RelSimplexDiameter(
    const std::vector<std::vector<T>>& S,
    const std::vector<int>& s,
    const T& epsilon
) {
    int d = s.size();
    T maxsofar = T{0};
    T temp = T{0};
    std::vector<std::vector<int>> subsets;
    std::vector<int> subset_t;
    std::vector<int> s2;
    std::vector<T> m(S[0].size());

    if (d > 2) {
        if (s[d - 1] >= S.size()) {
            s2 = s;
            s2.erase(s2.end() - 1);

            subset(s2, d - 1, 2, 0, subset_t, subsets);

            for (const auto& sigma : subsets) {
                auto first = S[sigma[0]].begin();
                auto first2 = S[sigma[1]].begin();

                for (size_t it = 0; it < m.size(); it++) {
                    m[it] = ((*first++) + (*first2++)) / 2;
                }

                if (vectorDistance(m.begin(), m.end(), S[0].begin()) < T{0.5} * epsilon) {
                    temp = BoundaryProjectedDistance(
                        S[sigma[0]].begin(), S[sigma[0]].end(),
                        S[sigma[1]].begin(),
                        S[0].begin(), S[0].end(),
                        T{0.5} * epsilon
                    );
                } else {
                    temp = T{0.5} * vectorDistance(S[sigma[0]].begin(), S[sigma[0]].end(), S[sigma[1]].begin());
                }

                if (temp > maxsofar) {
                    maxsofar = temp;
                }
            }
        } else {
            subset(s, d, 2, 0, subset_t, subsets);

            for (const auto& sigma : subsets) {
                temp = T{0.5} * vectorDistance(S[sigma[0]].begin(), S[sigma[0]].end(), S[sigma[1]].begin());

                if (temp > maxsofar) {
                    maxsofar = temp;
                }
            }
        }
    } else if (d > 1) {
        if (s[1] < S.size()) {
            maxsofar = T{0.5} * vectorDistance(S[s[0]].begin(), S[s[0]].end(), S[s[1]].begin());
        } else {
            temp = T{0.5} * epsilon - vectorDistance(S[0].begin(), S[0].end(), S[s[0]].begin());
            if (temp > T{0}) {
                maxsofar = temp;
            } else {
                maxsofar = T{0};
            }
        }
    }
    return maxsofar;
}


template <typename T>
void AddCofaces(
    const std::vector<std::vector<T>>& S,
    const std::vector<std::vector<int>>& G,
    const int& k,
    const std::vector<int>& t,
    const std::vector<int>& N,
    std::vector<std::vector<std::pair<T, std::vector<int>>>>& V
) {
    int d = t.size();
    V[d - 1].push_back({SimplexDiameter(S, t), t});

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
            AddCofaces(S, G, k, s, M, V);
        }
    }
}

template <typename T>
void AddRelCofaces(
    const std::vector<std::vector<T>>& S,
    const std::vector<std::vector<int>>& G,
    const int& k,
    const std::vector<int>& t,
    const std::vector<int>& N,
    std::vector<std::vector<std::pair<T, std::vector<int>>>>& V,
    const T& epsilon
) {
    int d = t.size();

    T diam = RelSimplexDiameter(S, t, epsilon);

    if (diam < T{0.5} * epsilon) {
        V[d - 1].push_back({diam, t});
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
                    std::set_intersection(
                        N.begin(), N.end(),
                        G[f].begin(), G[f].end(),
                        std::back_inserter(M)
                    );
                    AddRelCofaces(S, G, k, s, M, V, epsilon);
                } else {
                    AddRelCofaces(S, G, k, s, M, V, epsilon);
                }
            }
        }
    }
    return;
}

template <typename T>
void AddBndCofaces(
    const std::vector<std::vector<T>>& S,
    const std::vector<std::vector<int>>& G,
    const int& k,
    const std::vector<int>& t,
    const std::vector<int>& N,
    std::vector<std::vector<std::pair<T, std::vector<int>>>>& V,
    const T& epsilon
) {
    int d = t.size();
    V[d - 1].push_back({BndSimplexDiameter(S, t, epsilon), t});

    if (d < k + 2) {
        for (const int& f : N) {
            if (std::find(t.begin(), t.end(), f) == t.end()) {
                std::vector<int> M;
                std::vector<int> s = t;
                s.push_back(f);
                std::sort(s.begin(), s.end());
                std::set_intersection(
                    N.begin(), N.end(),
                    G[f].begin(), G[f].end(),
                    std::back_inserter(M)
                );
                AddBndCofaces(S, G, k, s, M, V, epsilon);
            }
        }
    }
    return;
}


template <typename T>
Cplx<T> IncrementalVR(const std::vector<std::vector<T>>& S, const T& epsilon, const int& k, const std::list<std::vector<std::vector<int>>>& nbh) {
    Cplx<T> VR(k + 2);
    VR.neighborhood = nbh;

    int n = VR.neighborhood.front().size();
    std::vector<int> u(1);

    for (int i = 0; i < n; i++) {
        u[0] = i;
        AddCofaces(S, VR.neighborhood.front(), k, u, VR.neighborhood.front()[i], VR.simplices);
    }

    int count = 0;

    for (int i = k + 1; i >= 0; --i) {
        std::stable_sort(VR.simplices[i].begin(), VR.simplices[i].end(), [](const std::pair<T, std::vector<int>>& left, const std::pair<T, std::vector<int>>& right) {
            return left.first < right.first;
        });

        for (auto simp_it = VR.simplices[i].rbegin(); simp_it != VR.simplices[i].rend(); simp_it++) {
            VR.IDs.emplace(simp_it->second, count);
            count += 1;
        }
    }

    return VR;
}

template <typename T>
Cplx<T> RelIncrementalVR(
    const std::vector<std::vector<T>>& S,
    const int& k,
    const std::list<std::vector<std::vector<int>>>& nbh,
    const T& epsilon
) {
    Cplx<T> VR(k + 2);
    VR.neighborhood = nbh;

    int n = VR.neighborhood.front().size();
    std::vector<int> u(1);

    for (int i = 0; i < n; i++) {
        if ((VR.neighborhood).front()[i].size() > 0) {
            u[0] = (VR.neighborhood).front()[i][0];
            AddRelCofaces(S, (VR.neighborhood).front(), k, u, (VR.neighborhood).front()[u[0]], VR.simplices, epsilon);
        }
    }

    int count = 0;

    for (int i = k + 1; i >= 0; --i) {
        std::stable_sort(VR.simplices[i].begin(), VR.simplices[i].end(), [](const std::pair<T, std::vector<int>>& left,
                                                                            const std::pair<T, std::vector<int>>& right) {
            return left.first < right.first;
        });

        for (auto simp_it = VR.simplices[i].rbegin(); simp_it != VR.simplices[i].rend(); simp_it++) {
            VR.IDs[simp_it->second] = count;
            count += 1;
        }
    }

    return VR;
}

template <typename T>
Cplx<T> BndIncrementalVR(
    const std::vector<std::vector<T>>& S,
    const int& k,
    const std::list<std::vector<std::vector<int>>>& nbh,
    const T& epsilon
) {
    Cplx<T> VR(k + 2);
    VR.neighborhood = nbh;

    int n = VR.neighborhood.front().size();
    std::vector<int> u(1);

    for (int i = 0; i < n; i++) {
        if ((VR.neighborhood).front()[i].size() > 0) {
            u[0] = (VR.neighborhood).front()[i][0];
            AddBndCofaces(S, (VR.neighborhood).front(), k, u, (VR.neighborhood).front()[u[0]], VR.simplices, epsilon);
        }
    }

    int count = 0;

    for (int i = k + 1; i >= 0; --i) {
        std::stable_sort(VR.simplices[i].begin(), VR.simplices[i].end(), [](const std::pair<T, std::vector<int>>& left,
                                                                            const std::pair<T, std::vector<int>>& right) {
            return left.first < right.first;
        });

        for (auto simp_it = VR.simplices[i].rbegin(); simp_it != VR.simplices[i].rend(); simp_it++) {
            VR.IDs[simp_it->second] = count;
            count += 1;
        }
    }

    return VR;
}
