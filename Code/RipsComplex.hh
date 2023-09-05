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


using namespace std;

struct Cplx
{
    Cplx()
    : simplices(), neighborhood(), IDs() {}
    vector< vector< pair<long double, vector<int> > > > simplices;
    list<vector< vector<int>>> neighborhood;
    ska::unordered_map<vector<int>, int, VectorHasher > IDs;

    Cplx(int n) {simplices = vector< vector< pair<long double,vector<int> > > >(n);};
};

long double SimplexDiameter
(
    const vector<vector<long double>>& S,
    const vector<int>& s
)
{
    int d = s.size();
    long double maxsofar = 0.0;
    long double temp = 0.0;

    if (d > 1)
    {
        for (int i = 0; i < d; i++)
        {
            for (int j = i+1; j < d; j++)
            {
                temp = vectorDistance(S[s[i]].begin(), S[s[i]].end(), S[s[j]].begin());
                
                if (temp > maxsofar)
                {
                    maxsofar = temp;
                }
            }
        }
    }
    return maxsofar;
}


long double BndSimplexDiameter
(
    const vector<vector<long double>>& S,
    const vector<int>& s,
    const long double& epsilon
)
{
    int d = s.size();
    long double maxsofar = 0.0;
    long double temp = 0.0;
    std::vector<long double> m(S[0].size());
    std::vector<long double>::const_iterator first;
    std::vector<long double>::const_iterator last;
    std::vector<long double>::const_iterator first2;

    if (d > 1)
    {
        for (int i = 0; i < d; i++)
        {
            for (int j = i+1; j < d; j++)
            {
                first = S[s[i]].begin();
                last = S[s[i]].end();
                first2 = S[s[j]].begin();

                for (int it = 0; it < m.size(); it++)
                {
                    m[it] = ((*first++) + (*first2++))/2;
                }

                if (vectorDistance(m.begin(),m.end(),S[0].begin()) < 0.5*epsilon)
                {
                    temp = G(S[s[i]],S[s[j]],S[0],0.5*epsilon);
                }

                else
                {
                    temp = 0.5*vectorDistance(S[s[i]].begin(), S[s[i]].end(), S[s[j]].begin());
                }
                   
                if (temp > maxsofar)
                {
                    maxsofar = temp;
                }
            }
        }
    }
    else
    {
        temp = 0.5*epsilon - vectorDistance(S[s[0]].begin(),S[s[0]].end(),S[0].begin());
        if (temp > 0)
        {
            maxsofar = temp;
        }
    }
    return maxsofar;
}


long double RelSimplexDiameter
(
    const vector<vector<long double>>& S,
    const vector<int>& s,
    const long double& epsilon
)
{
    int d = s.size();
    long double maxsofar = 0.0;
    long double temp = 0.0;
    vector<vector<int>> T;
    vector<int> t;
    vector<int> s2;
    std::vector<long double> m(S[0].size());
    std::vector<long double>::const_iterator first;
    std::vector<long double>::const_iterator last;
    std::vector<long double>::const_iterator first2;

    if (d > 2)
    {
        if (s[d-1] >= S.size())
        {
            s2 = s;

            s2.erase(s2.end()-1);

            subset(s2,d-1,2,0,t,T);
        
            for(const auto& sigma : T)
            {
                first = S[sigma[0]].begin();
                first2 = S[sigma[1]].begin();

                for (int it = 0; it < m.size(); it++)
                {
                    m[it] = ((*first++) + (*first2++))/2;
                }

                if (vectorDistance(m.begin(),m.end(),S[0].begin()) < 0.5*epsilon)
                {
                    temp = G(S[sigma[0]],S[sigma[1]],S[0],0.5*epsilon);
                }

                else
                {
                    temp = 0.5*vectorDistance(S[sigma[0]].begin(), S[sigma[0]].end(), S[sigma[1]].begin());
                }
                
                if (temp > maxsofar)
                {
                    maxsofar = temp; 
                }
            }
        }
        else
        {
            subset(s,d,2,0,t,T);
        
            for(const auto& sigma : T)
            {
                temp = 0.5*vectorDistance(S[sigma[0]].begin(),S[sigma[0]].end(),S[sigma[1]].begin());

                if (temp > maxsofar)
                {
                    maxsofar = temp; 
                }
            }
        }
    }
        
    else if(d > 1)
    {
        if (s[1] < S.size())
        {
            maxsofar = 0.5*vectorDistance(S[s[0]].begin(), S[s[0]].end(), S[s[1]].begin());
        }
        else
        {
            temp = (0.5*epsilon - vectorDistance(S[0].begin(), S[0].end(), S[s[0]].begin()));
            if (temp > 0)
                {
                    maxsofar = temp; 
                }
            else
            {
                maxsofar = 0.0;
            }
        }
    }
    return maxsofar;
}

void AddCofaces
(
const vector<vector<long double>>& S, 
const vector< vector<int>>& G,
const int& k,
const vector<int>& t,
const vector<int>& N,
vector<vector< pair<long double,vector<int> > > >& V
)
{
    int d = t.size();
    V[d-1].push_back({SimplexDiameter(S,t),t});

    vector<int> s;
    
    if (d > k+1)
    {
        return;
    }

    for(const int& f : N)
    {
        if ( (find(t.begin(), t.end(), f) != t.end()) == false )
        {
            vector<int> M;
            s = t;
            s.push_back(f);
            set_intersection(
            N.begin(), N.end(),
            G[f].begin(), G[f].end(),
            back_inserter( M ));
            AddCofaces(S, G, k, s, M, V);
        }
    }
}

void AddRelCofaces
(
const vector<vector<long double>>& S,
const vector< vector<int>>& G,
const int& k,
const vector<int>& t,
const vector<int>& N,
vector<vector< pair<long double,vector<int> > > >& V,
const long double epsilon
)
{
    int d = t.size();

    long double diam = RelSimplexDiameter(S,t,epsilon);

    if (diam < 0.5*epsilon)
    {
        V[d-1].push_back({diam,t});
    }
    else
    {
        return;
    }
    
    if (d < k+2)
    {
        for(const int& f : N)
        {
            if ( (find(t.begin(), t.end(), f) != t.end()) == false )
            {
                vector<int> M;
                vector<int> s = t;
                s.push_back(f);
                sort(s.begin(),s.end());

                if(s.back() < S.size())
                {
                    set_intersection(
                    N.begin(), N.end(),
                    G[f].begin(), G[f].end(),
                    back_inserter( M ));
                    AddRelCofaces(S, G, k, s, M, V, epsilon);
                }
                else
                {
                    AddRelCofaces(S, G, k, s, M, V, epsilon);
                }
                
            }
        }
    }
    return;
}

void AddBndCofaces
(
const vector<vector<long double>>& S,
const vector< vector<int>>& G,
const int& k,
const vector<int>& t,
const vector<int>& N,
vector<vector< pair<long double,vector<int> > > >& V,
const long double epsilon
)
{
    int d = t.size();
    V[d-1].push_back({BndSimplexDiameter(S,t,epsilon),t});

    if (d < k+2)
    {
        for(const int& f : N)
        {
            if ( (find(t.begin(), t.end(), f) != t.end()) == false )
            {
                vector<int> M;
                vector<int> s = t;
                s.push_back(f);
                sort(s.begin(),s.end());
                set_intersection(
                N.begin(), N.end(),
                G[f].begin(), G[f].end(),
                back_inserter( M ));
                AddBndCofaces(S, G, k, s, M, V, epsilon);   
            }
        }
    }
    return;
}

Cplx IncrementalVR
(
  const vector<vector<long double>>& S,
  const long double& epsilon,
  const int& k,
  const list<vector< vector<int>>>& nbh
)
{    
    Cplx VR(k+2);
    VR.neighborhood = nbh;

    int n = VR.neighborhood.front().size();

    vector<int> u(1);

    for(int i = 0; i<n; i++)
    {
        u[0] = i;
        AddCofaces(S, (VR.neighborhood).front(), k, u, (VR.neighborhood).front()[i], VR.simplices);
    }

    int count = 0;

    for (int i = k+1; i>= 0; --i)
    {
        stable_sort(VR.simplices[i].begin(), VR.simplices[i].end(), [](const pair<long double,vector<int> >&  left,
        const pair<long double,vector<int> >&  right) {
        return left.first < right.first;}
    );
        for(vector<pair<long double,vector<int>>>::reverse_iterator simp_it = (VR.simplices[i]).rbegin(); simp_it != VR.simplices[i].rend(); simp_it++)
        {
            VR.IDs.emplace(simp_it->second,count);
            count += 1;
        }
    }

return VR;
}

Cplx RelIncrementalVR
(
  const vector<vector<long double>>& S,
  const int& k,
  const list<vector< vector<int>>>& nbh,
  const long double epsilon
)
{    
    Cplx VR(k+2);
    VR.neighborhood = nbh;

    int n = VR.neighborhood.front().size();

    vector<int> u(1);

    for(int i = 0; i<n; i++)
    {
        if((VR.neighborhood).front()[i].size() > 0)
        {
            u[0] = (VR.neighborhood).front()[i][0];
            AddRelCofaces(S, (VR.neighborhood).front(), k, u, (VR.neighborhood).front()[u[0]], VR.simplices, epsilon);
        }
    }

    int count = 0;

    for (int i = k+1; i>= 0; --i)
    {
        stable_sort(VR.simplices[i].begin(), VR.simplices[i].end(), [](const pair<long double,vector<int> >&  left,
        const pair<long double,vector<int> >&  right) {
        return left.first < right.first;});

        for(vector<pair<long double,vector<int>>>::reverse_iterator simp_it = (VR.simplices[i]).rbegin(); simp_it != VR.simplices[i].rend(); simp_it++)
        {
            VR.IDs[simp_it->second] = count;
            count += 1;
        }
    }

return VR;
}

Cplx BndIncrementalVR
(
  const vector<vector<long double>>& S,
  const int& k,
  const list<vector< vector<int>>>& nbh,
  const long double epsilon
)
{    
    Cplx VR(k+2);
    VR.neighborhood = nbh;

    int n = VR.neighborhood.front().size();

    vector<int> u(1);

    for(int i = 0; i<n; i++)
    {
        if((VR.neighborhood).front()[i].size() > 0)
        {
            u[0] = (VR.neighborhood).front()[i][0];
            AddBndCofaces(S, (VR.neighborhood).front(), k, u, (VR.neighborhood).front()[u[0]], VR.simplices, epsilon);
        }
    }

    int count = 0;

    for (int i = k+1; i>= 0; --i)
    {
        stable_sort(VR.simplices[i].begin(), VR.simplices[i].end(), [](const pair<long double,vector<int> >&  left,
        const pair<long double,vector<int> >&  right) {
        return left.first < right.first;});

        for(vector<pair<long double,vector<int>>>::reverse_iterator simp_it = (VR.simplices[i]).rbegin(); simp_it != VR.simplices[i].rend(); simp_it++)
        {
            VR.IDs[simp_it->second] = count;
            count += 1;
        }
    }

return VR;
}
