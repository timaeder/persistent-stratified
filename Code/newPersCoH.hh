#pragma once
#include <iostream>
#include <fstream>
#include <list>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <queue>
#include <chrono>
#include <iomanip>
#include "newRipsComplex.hh"
#include "newNeighborhood.hh"

using namespace std;

// Adjusted SimpCoBound function
vector<int> SimpCoBound(
    const vector<int>& s,
    const vector<vector<int>>& N,
    const ska::unordered_map<vector<int>, int, VectorHasher>& IDs,
    const int& n
) {
    int k = s.size() - 1;
    vector<int> cobound;
    vector<int> last_intersection = N[s[0]];
    vector<int> curr_intersection;

    if (last_intersection.size() > k + 1) {
        bool stop = false;
        for (int i = 1; i < s.size(); i++) {
            if (N[s[i]].size() > 0 && !stop) {
                set_intersection(
                    last_intersection.begin(), last_intersection.end(),
                    N[s[i]].begin(), N[s[i]].end(),
                    back_inserter(curr_intersection)
                );
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
                    if (find(s.begin(), s.end(), v) == s.end()) {
                        vector<int> t = s;
                        t.push_back(v);
                        sort(t.begin(), t.end());
                        cobound.push_back(IDs.find(t)->second);
                    }
                }
            }
        }
        if (cobound.size()) {
            sort(cobound.begin(), cobound.end());
        }
    }
    return cobound;
}

// Adjusted RelSimpCoBound function
vector<int> RelSimpCoBound(
    const vector<vector<long double>>& S,
    const vector<int>& s,
    const vector<vector<int>>& N,
    const ska::unordered_map<vector<int>, int, VectorHasher>& IDs,
    const int& n,
    const long double epsilon
) {
    int k = s.size() - 1;
    vector<int> cobound;
    vector<int> last_intersection = N[s[0]];
    vector<int> curr_intersection;
    int conept = S.size();

    if (last_intersection.size() > k + 1) {
        bool stop = false;
        for (int i = 1; i < s.size(); i++) {
            if (N[s[i]].size() > 0 && !stop) {
                set_intersection(
                    last_intersection.begin(), last_intersection.end(),
                    N[s[i]].begin(), N[s[i]].end(),
                    back_inserter(curr_intersection)
                );
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
                    if (find(s.begin(), s.end(), v) == s.end()) {
                        vector<int> t = s;
                        t.push_back(v);
                        sort(t.begin(), t.end());
                        if (t.back() != conept) {
                            cobound.push_back(IDs.find(t)->second);
                        } else {
                            Complex complex(n + 1);
                            if (complex.computeRelSimplexDiameter(S, t, epsilon) < 0.5 * epsilon) {
                                cobound.push_back(IDs.find(t)->second);
                            }
                        }
                    }
                }
            }
        }
        if (cobound.size()) {
            sort(cobound.begin(), cobound.end());
        }
    }
    return cobound;
}

// Adjusted AddBndry function
void AddBndry(vector<int>& cbID1, vector<int>& cbID2) {
    vector<int> symm_diff;

    set_symmetric_difference(
        cbID1.begin(), cbID1.end(),
        cbID2.begin(), cbID2.end(),
        back_inserter(symm_diff)
    );

    std::swap(cbID1, symm_diff);
}

// Adjusted Copairings function
vector<vector<pair<long double, long double>>> Copairings(
    const vector<vector<long double>>& S,
    Complex& VR,
    const long double& epsilon,
    const int& k
) {
    int n = VR.simplices.size(); // Access simplices directly
    vector<vector<pair<long double, long double>>> copairs(n);

    int full_size = VR.IDs.size();
    int curr_size = 0;

    vector<vector<int>> redcols(full_size, vector<int>{-1});
    vector<set<int>> redrows(full_size);

    for (int d = 0; d < n - 1; d++) {
        curr_size += VR.simplices[d].size(); // Access simplices directly

        auto simp_it1 = VR.simplices[d].rbegin(); // Reverse iterator for dimension d
        while (simp_it1 != VR.simplices[d].rend()) {
            const vector<int> simp1 = simp_it1->second;
            int id1 = full_size - curr_size + (simp_it1 - VR.simplices[d].rbegin());

            if (redcols[id1].size()) {
                auto cbID1 = (redcols[id1] = SimpCoBound(simp1, VR.neighborhood.back(), VR.IDs, n - 1));
                if (cbID1.size()) {
                    int pivot1 = cbID1.back();
                    bool stop = false;

                    while (redrows[pivot1].size() > 0 && !stop) {
                        if (redrows[pivot1].size()) {
                            auto cbID2 = redcols[*(redrows[pivot1].begin())];
                            AddBndry(cbID1, cbID2);

                            if (cbID1.size()) {
                                pivot1 = cbID1.back();
                            } else {
                                stop = true;
                            }
                        } else {
                            stop = true;
                        }
                    }
                }

                redcols[id1] = cbID1;

                if (cbID1.size()) {
                    redrows[cbID1.back()].insert(id1);
                    if (redcols[*(cbID1.rbegin())].size()) {
                        long double pers = std::fabs(
                            VR.simplices[d + 1][full_size - curr_size - *(cbID1.rbegin()) - 1].first
                            - simp_it1->first
                        );

                        if (pers > 0.001) {
                            copairs[d].push_back(make_pair(simp_it1->first, pers));
                        }
                        redcols[*(cbID1.rbegin())].clear();
                    }
                } else {
                    copairs[d].push_back(make_pair(simp_it1->first, epsilon + 0.1));
                }
            }
            ++simp_it1;
        }
    }

    return copairs;
}

vector<vector<pair<long double, long double>>> CLocGCopairings(
    const vector<vector<long double>>& S,
    const int& k,
    const vector<long double>& z,
    const long double rad
) {
    vector<pair<long double, vector<long double>>> FiltPts;
    long double epsilon = rad;
    long double r;

    // Filter points within the radius
    for (const auto& v : S) {
        r = DistanceCalculator::vectorDistance(z.begin(), z.end(), v.begin());
        if (r <= rad) {
            FiltPts.push_back(make_pair(r, v));
        }
    }

    // Sort filtered points by distance
    stable_sort(FiltPts.begin(), FiltPts.end());

    vector<vector<long double>> LPC;
    for (const auto& v : FiltPts) {
        LPC.push_back(v.second);
    }

    // Create Complex and Neighborhood objects
    Complex VR(k + 2);
    Neighborhood neighborhoodCalculator(LPC, rad);

    // Compute CVR neighborhood and RelIncrementalVR
    auto nbh = neighborhoodCalculator.computeCVRNeighborhood();
    VR.computeRelIncrementalVR(LPC, rad, k);

    int n = k + 1;
    for (int i = 0; i <= k + 1; i++) {
        if (!VR.simplices[i].size()) {
            n = i - 1;
            break;
        }
    }

    int full_size = VR.IDs.size();
    int curr_size = 0;

    vector<vector<int>> redcols(full_size, vector<int>{-1});
    vector<set<int>> redrows(full_size);

    vector<vector<pair<long double, long double>>> LocCopairs(k + 1);

    for (int d = 0; d < n; d++) {
        curr_size += VR.simplices[d].size();

        auto simp_it1 = VR.simplices[d].rbegin();
        while (simp_it1 != VR.simplices[d].rend()) {
            const vector<int> simp1 = simp_it1->second;
            int id1 = full_size - curr_size + (simp_it1 - VR.simplices[d].rbegin());

            if (redcols[id1].size()) {
                auto cbID1 = (redcols[id1] = RelSimpCoBound(LPC, simp1, VR.neighborhood.back(), VR.IDs, n, rad));
                if (cbID1.size()) {
                    int pivot1 = cbID1.back();
                    bool stop = false;

                    while (redrows[pivot1].size() > 0 && !stop) {
                        if (redrows[pivot1].size()) {
                            auto cbID2 = redcols[*(redrows[pivot1].begin())];
                            AddBndry(cbID1, cbID2);

                            if (cbID1.size()) {
                                pivot1 = cbID1.back();
                            } else {
                                stop = true;
                            }
                        } else {
                            stop = true;
                        }
                    }
                }

                redcols[id1] = cbID1;

                if (cbID1.size()) {
                    redrows[cbID1.back()].insert(id1);
                    if (redcols[*(cbID1.rbegin())].size()) {
                        long double pers = std::fabs(
                            VR.simplices[d + 1][full_size - curr_size - *(cbID1.rbegin()) - 1].first
                            - simp_it1->first
                        );

                        if (pers > 0.001) {
                            LocCopairs[d].push_back(make_pair(simp_it1->first, pers));
                        }
                        redcols[*(cbID1.rbegin())].clear();
                    }
                } else {
                    LocCopairs[d].push_back(make_pair(simp_it1->first, (0.5 * rad) - simp_it1->first + 0.001));
                }
            }
            ++simp_it1;
        }
    }

    return LocCopairs;
}

// Adjusted BndLocGCopairings function
vector<vector<pair<long double, long double>>> BndLocGCopairings(
    const vector<vector<long double>>& S,
    const int& k,
    const vector<long double>& z,
    const long double rad
) {
    vector<pair<long double, vector<long double>>> FiltPts;
    long double epsilon = rad;
    long double r;

    // Filter points within the radius
    for (const auto& v : S) {
        r = DistanceCalculator::vectorDistance(z.begin(), z.end(), v.begin());
        if (r <= rad) {
            FiltPts.push_back(make_pair(r, v));
        }
    }

    // Sort filtered points by distance
    stable_sort(FiltPts.begin(), FiltPts.end());

    vector<vector<long double>> LPC;
    for (const auto& v : FiltPts) {
        LPC.push_back(v.second);
    }

    // Create Complex and Neighborhood objects
    Complex VR(k + 2);
    Neighborhood neighborhoodCalculator(LPC, rad);

    // Compute boundary neighborhood and BndIncrementalVR
    auto nbh = neighborhoodCalculator.computeBndNeighborhood();
    VR.computeBndIncrementalVR(LPC, rad, k);

    int n = k + 1;
    for (int i = 0; i <= k + 1; i++) {
        if (!VR.simplices[i].size()) {
            n = i - 1;
            break;
        }
    }

    int full_size = VR.IDs.size();
    int curr_size = 0;

    vector<vector<int>> redcols(full_size, vector<int>{-1});
    vector<set<int>> redrows(full_size);

    vector<vector<pair<long double, long double>>> LocCopairs(k + 1);

    for (int d = 0; d < n; d++) {
        curr_size += VR.simplices[d].size();

        auto simp_it1 = VR.simplices[d].rbegin();
        while (simp_it1 != VR.simplices[d].rend()) {
            const vector<int> simp1 = simp_it1->second;
            int id1 = full_size - curr_size + (simp_it1 - VR.simplices[d].rbegin());

            if (redcols[id1].size()) {
                auto cbID1 = (redcols[id1] = SimpCoBound(simp1, VR.neighborhood.back(), VR.IDs, n));
                if (cbID1.size()) {
                    int pivot1 = cbID1.back();
                    bool stop = false;

                    while (redrows[pivot1].size() > 0 && !stop) {
                        if (redrows[pivot1].size()) {
                            auto cbID2 = redcols[*(redrows[pivot1].begin())];
                            AddBndry(cbID1, cbID2);

                            if (cbID1.size()) {
                                pivot1 = cbID1.back();
                            } else {
                                stop = true;
                            }
                        } else {
                            stop = true;
                        }
                    }
                }

                redcols[id1] = cbID1;

                if (cbID1.size()) {
                    redrows[cbID1.back()].insert(id1);
                    if (redcols[*(cbID1.rbegin())].size()) {
                        long double pers = std::fabs(
                            VR.simplices[d + 1][full_size - curr_size - *(cbID1.rbegin()) - 1].first
                            - simp_it1->first
                        );

                        if (pers > 0.001) {
                            LocCopairs[d].push_back(make_pair(simp_it1->first, pers));
                        }
                        redcols[*(cbID1.rbegin())].clear();
                    }
                } else {
                    LocCopairs[d].push_back(make_pair(simp_it1->first, (0.5 * rad) - simp_it1->first + 0.001));
                }
            }
            ++simp_it1;
        }
    }

    // Reduce homology
    LocCopairs[0].erase(LocCopairs[0].end() - 1);

    return LocCopairs;
}