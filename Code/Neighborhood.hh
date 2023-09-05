#pragma once
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <algorithm>
#include "EuclDist.hh"

std::list<std::vector< std::vector<int>>> Neighborhood
( 
    const std::vector<std::vector<long double>>& S,
    const long double& epsilon
)
{
    int n = S.size();
    long double dist;
    std::vector<std::vector<int>> graph1(n);
    std::vector<std::vector<int>> graph2(n);

    for(int i = 0; i < n; i++)
    {
        int j = 0;
        while (j < n)
        {
            dist = vectorDistance(S[i].begin(),S[i].end(),S[j].begin());
            if (dist <= epsilon)
            {
                if (j >= i)
                {
                    graph1[i].push_back(j);
                }
                graph2[i].push_back(j);
                
            }
            j = j+1;
        }
    }
    std::list<std::vector< std::vector<int>>> graphs;
    graphs.push_back(graph1);
    graphs.push_back(graph2);
    return graphs;
} 

std::list<std::vector< std::vector<int>>> CVRNeighborhood
( 
    const std::vector<std::vector<long double>>& LS,
    const long double& epsilon
)
{
    int k = LS.size();
    long double dist = 0.0;
    std::vector<std::set<int>> graph1(k+1);
    std::vector<std::set<int>> graph2(k+1);
    
    for(int i = 0; i < k; i++)
    {
        graph2[i].insert(i);
        graph1[i].insert(i);

        for (int j = i+1; j < k; j++)
        {
            dist = 0.5*vectorDistance(LS[i].begin(),LS[i].end(),LS[j].begin());
            if (dist < 0.5*epsilon)
            {
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

    for (const auto& g : graph1)
    {   
        graph3.push_back(std::vector<int> (g.begin(),g.end()));
    }

    for (const auto& g : graph2)
    {
        graph4.push_back(std::vector<int> (g.begin(),g.end()));
    }
    std::list<std::vector< std::vector<int>>> graphs;
    graphs.push_back(graph3);
    graphs.push_back(graph4);
    return graphs;
} 

std::list<std::vector< std::vector<int>>> BndNeighborhood
( 
    const std::vector<std::vector<long double>>& LS,
    const long double& epsilon
)
{
    int k = LS.size();
    long double dist = 0.0;
    std::vector<long double> m(LS[0].size());
    std::vector<long double>::const_iterator first;
    std::vector<long double>::const_iterator last;
    std::vector<long double>::const_iterator first2;
    std::vector<std::set<int>> graph1(k);
    std::vector<std::set<int>> graph2(k);

    for(int i = 0; i < k; i++)
    {
        if ((0.5*epsilon - vectorDistance(LS[i].begin(),LS[i].end(),LS[0].begin())) < 0.5*epsilon)
        {

        graph2[i].insert(i);
        graph1[i].insert(i);
        for (int j = i+1; j < k; j++)
        {
            first = LS[i].begin();
            last = LS[i].end();
            first2 = LS[j].begin();

            for (int it = 0; it < m.size(); it++)
            {
                m[it] = ((*first++) + (*first2++))/2;
            }

            if (vectorDistance(m.begin(),m.end(),LS[0].begin()) < 0.5*epsilon)
            {
                dist = G(LS[i],LS[j],LS[0],0.5*epsilon);
            }

            else
            {
                dist = 0.5*vectorDistance(LS[i].begin(),LS[i].end(),LS[j].begin());
            }

            if (dist < 0.5*epsilon)
            {
                graph1[i].insert(j);
                graph2[i].insert(j);
                graph2[j].insert(i);
            }
        }
        }
    }

    std::vector<std::vector<int>> graph3;
    std::vector<std::vector<int>> graph4;

    for (const auto& g : graph1)
    {
        graph3.push_back(std::vector<int> (g.begin(),g.end()));
    }

    for (const auto& g : graph2)
    {
        graph4.push_back(std::vector<int> (g.begin(),g.end()));
    }
    std::list<std::vector< std::vector<int>>> graphs;
    graphs.push_back(graph3);
    graphs.push_back(graph4);
    return graphs;
} 
