Hi I need your help with modifying some existing c++ code I have. I started transforming certiain scripts to a more object oriented iteration but one script is not adjusted to that and still uses the older version of the already transformed code. I would appreciate if you could take a look at the already transformed code I will show your first and then look at the code to be changed and make the adjustements so that it is also object oriented and works with the other code. First here is my already updated code here the class of simplicial complex is the class that should be extended by new methods from the code to be changed: #pragma once
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
};#pragma once
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
}; And here is the code that needs to be changed and implemented as a class method of simplicial complex class: #pragma once
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
#include "RipsComplex.hh"

using namespace std;

vector< int > SimpCoBound
( 
    const vector<int>& s, 
    const vector< vector<int>>& N,
    const ska::unordered_map<vector<int>, int, VectorHasher>& IDs,
    const int& n
)
{
    int k = s.size()-1;
    vector< int > cobound;
    vector<int> last_intersection = N[s[0]];
    vector<int> curr_intersection;

    if (last_intersection.size()>k+1)
    {
        bool stop = false;
        for (int i = 1; i < s.size(); i++)
        {
            if (N[s[i]].size() > 0 && !(stop))
            {
                set_intersection(last_intersection.begin(), last_intersection.end(),
                    (N[s[i]]).begin(), (N[s[i]]).end(),
                    back_inserter(curr_intersection));
                std::swap(last_intersection, curr_intersection);
                curr_intersection.clear();
            }
            else
            {
                stop = true;
                last_intersection.clear();
            }
            
        }
        if (k < n)
        {
            if(last_intersection.size())
            {
                for (const int& v : last_intersection)
                {
                    if ((find(s.begin(), s.end(), v) != s.end()) == false )
                    {
                        vector<int> t;
                        t = s;
                        t.push_back(v);
                        sort(t.begin(),t.end());
                        cobound.push_back(IDs.find(t)->second);
                    }
                }
            }
        }
        if(cobound.size())
        {sort(cobound.begin(),cobound.end());}
    }
    return cobound;
}

vector< int > RelSimpCoBound
( 
    const vector<vector<long double>>& S,
    const vector<int>& s, 
    const vector< vector<int>>& N,
    const ska::unordered_map<vector<int>, int, VectorHasher>& IDs,
    const int& n,
    const long double epsilon
)
{

    int k = s.size()-1;
    vector<int> cobound;
    vector<int> last_intersection = N[s[0]];
    vector<int> curr_intersection;
    int conept = S.size();

    if (last_intersection.size()>k+1)
    {
        bool stop = false;
        for (int i = 1; i < s.size(); i++)
        {
            if (N[s[i]].size() > 0 && !(stop))
            {
                set_intersection(last_intersection.begin(), last_intersection.end(),
                    (N[s[i]]).begin(), (N[s[i]]).end(),
                    back_inserter(curr_intersection));
                std::swap(last_intersection, curr_intersection);
                curr_intersection.clear();
            }
            else
            {
                stop = true;
                last_intersection.clear();
            }
            
        }
        if (k < n)
        {
            if(last_intersection.size())
            {
                for (const int& v : last_intersection)
                {
                    if ( (find(s.begin(),
                         s.end(), v) != s.end()) == false )
                    {
                        vector<int> t = s;
                        t.push_back(v);
                        sort(t.begin(),t.end());
                        if (t.back() != conept)
                        {
                            cobound.push_back(IDs.find(t)->second);
                        }
                        else
                        {
                            if(RelSimplexDiameter(S,t,epsilon) < 0.5*epsilon)
                            {cobound.push_back(IDs.find(t)->second);}
                        }   
                    }
                }
            }
        }
        if(cobound.size())
        {sort(cobound.begin(),cobound.end());}
    }
    return cobound;
}

void AddBndry(vector<int>& cbID1, vector<int>& cbID2)
{
    vector<int> symm_diff;
    
    set_symmetric_difference
    (
        cbID1.begin(), cbID1.end(),
        cbID2.begin(), cbID2.end(),
        back_inserter(symm_diff)
    );
    
    std::swap(cbID1,symm_diff);
}

vector < vector< pair< long double, long double > > > Copairings
(
    const vector<vector<long double>>& S,
    const Cplx& VR,
    const long double& epsilon,
    const int& k, 
    const list<vector< vector<int>>>& nbh
)
{
    int n = VR.simplices.size();
    int count;
    int step = 0;
    bool stop;
    bool test;

    vector < vector< pair< long double, long double > > > copairs(n);

    int id1;
    int id2;
    vector<int> cbID1;
    vector<int> cbID2;
    int pivot1;

    int full_size = VR.IDs.size();
    int curr_size = 0;
    
    
    vector< vector<int> > redcols(full_size);
    for (int j = 0; j < full_size; j++)
    {
        redcols[j] = vector<int>{-1};
    }

    for (int d = 0; d < n-1; d++)
    {
        vector< set<int> > redrows(full_size);
        curr_size = curr_size + VR.simplices[d].size();

        std::vector< pair<long double,vector<int> > >::const_reverse_iterator simp_it1;    
        simp_it1 = (VR.simplices[d]).rbegin();

        while (simp_it1 != VR.simplices[d].rend())
        {
            const vector<int> simp1 = simp_it1->second;
            id1 = full_size - curr_size + (simp_it1-VR.simplices[d].rbegin());
            
            if ( redcols[id1].size() )
            {
                cbID1 = (redcols[id1] = SimpCoBound(simp1,VR.neighborhood.back(),VR.IDs,n-1));
                if(cbID1.size())
                {
                    pivot1 = cbID1.back();

                    stop = false;

                    while( redrows[pivot1].size() > 0 && !(stop) )
                    {
                        if (redrows[pivot1].size())
                        {
                            cbID2 = redcols[*(redrows[pivot1].begin())];

                            AddBndry(cbID1,cbID2);

                            if (cbID1.size())
                            {
                                pivot1 = cbID1.back();
                            }

                            else
                            {
                                stop = true;
                            }
                        }
                        else
                        {
                            stop = true;
                        }
                            
                    }
                }
                
                redcols[id1] = cbID1;

                if (cbID1.size())
                {
                    redrows[pivot1].insert(id1);
                    if( redcols[*(cbID1.rbegin())].size() )
                    {
                        if ( std::abs(VR.simplices[d+1][full_size - curr_size - *(cbID1.rbegin())-1].first - simp_it1->first) > 0.0001 )
                        {
                            copairs[d].push_back(make_pair(simp_it1->first,VR.simplices[d+1][full_size - curr_size - *(cbID1.rbegin())-1].first));
                        }
                        redcols[*(cbID1.rbegin())].clear();
                    }
                }
                else
                {
                    copairs[d].push_back(make_pair(simp_it1->first,epsilon+0.1));
                }
            }
            ++simp_it1;
        }
    }

    return copairs;
}


vector < vector< pair<long double, long double > > > CLocGCopairings
(
    const vector<vector<long double>>& S,
    const int& k,
    const vector<long double>& z,
    const long double rad
)
{

    vector< pair<long double, vector<long double> > > FiltPts;
    long double epsilon = rad;
    long double r;

    for (const auto& v : S)
    {
        r = vectorDistance(z.begin(),z.end(),v.begin());
        if (r <= rad)
        {
            FiltPts.push_back(make_pair(r,v));
        }
    }

    stable_sort(FiltPts.begin(),FiltPts.end());

    vector<vector<long double>> LPC;

    for (const auto& v : FiltPts)
    {
        LPC.push_back(v.second);
    }

    Cplx VR;

    list<vector< vector<int>>> nbh;

    int count;
    int step;
    bool stop;
    bool test;
    long double pers = 0.0;

    int id1;
    int id2;
    vector<int> cbID1;
    vector<int> cbID2;
    int pivot1;

    int full_size;
    int curr_size;

    vector< vector<int> > redcols(full_size);
    vector< set<int> > redrows(full_size);

    vector< vector < pair<long double, long double > > > LocCopairs(k+1);

    nbh = CVRNeighborhood(LPC, rad);

    VR = RelIncrementalVR(LPC, k, nbh, rad);

    int n = k+1;

    for (int i = 0; i <= k+1; i ++)
    {
        if( !VR.simplices[i].size() )
        {
            n = i-1;
            break;
        }               
    }

    step = 0;

    full_size = VR.IDs.size();
    curr_size = 0;
        
    for (int j = 0; j < full_size; j++)
    {
        redcols.push_back(vector<int>{-1});
        redrows.push_back(set<int>{});
    }


    for (int d = 0; d < n; d++)
    {
        curr_size = curr_size + VR.simplices[d].size();

        std::vector< pair<long double,vector<int> > >::const_reverse_iterator simp_it1;
        simp_it1 = (VR.simplices[d]).rbegin();

            while (simp_it1 != VR.simplices[d].rend())
            {

                const vector<int> simp1 = simp_it1->second;
                id1 = full_size - curr_size + (simp_it1-VR.simplices[d].rbegin());

                if (redcols[id1].size())
                {
                    cbID1 = (redcols[id1] = RelSimpCoBound(LPC,simp1,VR.neighborhood.back(),VR.IDs,n,rad));
                    if (cbID1.size())
                    {
                        pivot1 = cbID1.back();

                        stop = false;

                        while( redrows[pivot1].size() > 0 && !(stop) )
                        {
                            if (redrows[pivot1].size())
                            {
                                cbID2 = redcols[*(redrows[pivot1].begin())];

                                AddBndry(cbID1,cbID2);           

                                if (cbID1.size())
                                {
                                    pivot1 = cbID1.back();
                                }

                                else
                                {
                                    stop = true;
                                }
                            }
                            else
                            {
                                stop = true;
                            }
                                    
                        }
                    }
                                            
                    redcols[id1] = cbID1;

                    if (cbID1.size())
                    {
                        redrows[pivot1].insert(id1);
                        if( redcols[*(cbID1.rbegin())].size())
                        {
                            pers = std::fabs(
                                VR.simplices[d+1][full_size - curr_size - *(cbID1.rbegin()) -1].first
                                - simp_it1->first);

                            if (pers > 0.001)
                            {
                                LocCopairs[d].push_back(
                                make_pair(simp_it1->first, pers));
                            }
                            redcols[*(cbID1.rbegin())].clear();
                        }
                    }

                    else
                    {
                        LocCopairs[d].push_back(make_pair(simp_it1->first,(0.5*rad) - simp_it1->first + 0.001));
                    }
                }
                ++simp_it1;
            }
    }

    return LocCopairs;
}

vector < vector< pair<long double, long double > > > BndLocGCopairings
(
    const vector<vector<long double>>& S,
    const int& k,
    const vector<long double>& z,
    const long double rad
)
{

    vector< pair<long double, vector<long double> > > FiltPts;
    long double epsilon = rad;
    long double r;

    for (const auto& v : S)
    {
        r = vectorDistance(z.begin(),z.end(),v.begin());
        if (r <= rad)
        {
            FiltPts.push_back(make_pair(r,v));
        }
    }

    stable_sort(FiltPts.begin(),FiltPts.end());

    vector<vector<long double>> LPC;

    for (const auto& v : FiltPts)
    {
        LPC.push_back(v.second);
    }

    Cplx VR;

    list<vector< vector<int>>> nbh;

    int count;
    int step;
    bool stop;
    bool test;
    long double pers = 0.0;

    int id1;
    int id2;
    vector<int> cbID1;
    vector<int> cbID2;
    int pivot1;

    int full_size;
    int curr_size ;

    vector< vector<int> > redcols(full_size);
    vector< set<int> > redrows(full_size);

    vector< vector < pair<long double, long double > > > LocCopairs(k+1);

    nbh = BndNeighborhood(LPC, rad);

    VR = BndIncrementalVR(LPC, k, nbh, rad);

    int n = k+1;

    for (int i = 0; i <= k+1; i ++)
    {
        if( !VR.simplices[i].size() )
        {
            n = i-1;
            break;
        }               
    }

    step = 0;

    full_size = VR.IDs.size();
    curr_size = 0;
        
    for (int j = 0; j < full_size; j++)
    {
        redcols.push_back(vector<int>{-1});
        redrows.push_back(set<int>{});
    }


    for (int d = 0; d < n; d++)
    {
        curr_size = curr_size + VR.simplices[d].size();

        std::vector< pair<long double,vector<int> > >::const_reverse_iterator simp_it1;
        simp_it1 = (VR.simplices[d]).rbegin();

        while (simp_it1 != VR.simplices[d].rend())
        {
            const vector<int> simp1 = simp_it1->second;
            id1 = full_size - curr_size + (simp_it1-VR.simplices[d].rbegin());
                
            if (redcols[id1].size())
            {
                cbID1 = (redcols[id1] = SimpCoBound(simp1,VR.neighborhood.back(),VR.IDs,n));
                if (cbID1.size())
                {
                    pivot1 = cbID1.back();

                    stop = false;

                    while( redrows[pivot1].size() > 0 && !(stop) )
                    {
                        if (redrows[pivot1].size())
                        {
                            cbID2 = redcols[*(redrows[pivot1].begin())];

                            AddBndry(cbID1,cbID2);           

                            if (cbID1.size())
                            {
                                pivot1 = cbID1.back();
                            }

                            else
                            {
                                stop = true;
                            }
                        }
                        else
                        {
                            stop = true;
                        }
                                
                    }
                }
                                        
                redcols[id1] = cbID1;

                if (cbID1.size())
                {
                    redrows[pivot1].insert(id1);
                    if( redcols[*(cbID1.rbegin())].size())
                    {
                        pers = std::abs(
                            VR.simplices[d+1][full_size - curr_size - *(cbID1.rbegin()) - 1].first
                            - simp_it1->first);

                        if (pers > 0.001)
                        {
                            LocCopairs[d].push_back(
                            make_pair(simp_it1->first, pers));
                        }
                        redcols[*(cbID1.rbegin())].clear();
                    }
                }

                else
                {
                    LocCopairs[d].push_back(make_pair(simp_it1->first,(0.5*rad) - simp_it1->first+0.001));
                }
            }
            ++simp_it1;
        }

    }

    //"reduce homology"
    LocCopairs[0].erase(LocCopairs[0].end()-1);

    return LocCopairs;
}