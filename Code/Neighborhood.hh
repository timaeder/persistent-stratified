#pragma once
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <cmath>
#include <algorithm>
#include "EuclDist.hh"

template <typename T>
std::list<std::vector<std::vector<int>>> Neighborhood(
    const std::vector<std::vector<T>>& S,
    const T& epsilon
)
{
    int n = S.size();
    T dist;
    std::vector<std::vector<int>> graph1(n);
    std::vector<std::vector<int>> graph2(n);

    for (int i = 0; i < n; i++)
    {
        int j = 0;
        while (j < n)
        {
            dist = vectorDistance(S[i].begin(), S[i].end(), S[j].begin());
            if (dist <= epsilon)
            {
                if (j >= i)
                {
                    graph1[i].push_back(j);
                }
                graph2[i].push_back(j);
            }
            j = j + 1;
        }
    }
    std::list<std::vector<std::vector<int>>> graphs;
    graphs.push_back(graph1);
    graphs.push_back(graph2);
    return graphs;
}

template <typename T>
std::list<std::vector<std::vector<int>>> CVRNeighborhood(
    const std::vector<std::vector<T>>& LS,
    const T& epsilon
)
{
    int k = LS.size();
    T dist = static_cast<T>(0.0);
    std::vector<std::set<int>> graph1(k + 1);
    std::vector<std::set<int>> graph2(k + 1);

    for (int i = 0; i < k; i++)
    {
        graph2[i].insert(i);
        graph1[i].insert(i);

        for (int j = i + 1; j < k; j++)
        {
            dist = static_cast<T>(0.5) * vectorDistance(LS[i].begin(), LS[i].end(), LS[j].begin());
            if (dist < static_cast<T>(0.5) * epsilon)
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

    for (typename std::vector<std::set<int>>::const_iterator it = graph1.begin(); it != graph1.end(); ++it)
    {
        graph3.push_back(std::vector<int>(it->begin(), it->end()));
    }

    for (typename std::vector<std::set<int>>::const_iterator it = graph2.begin(); it != graph2.end(); ++it)
    {
        graph4.push_back(std::vector<int>(it->begin(), it->end()));
    }

    std::list<std::vector<std::vector<int>>> graphs;
    graphs.push_back(graph3);
    graphs.push_back(graph4);
    return graphs;
}

template <typename T>
std::list<std::vector<std::vector<int>>> BndNeighborhood(
    const std::vector<std::vector<T>>& LS,
    const T& epsilon
)
{
    int k = LS.size();
    T dist = static_cast<T>(0.0);
    std::vector<T> m(LS[0].size());
    std::vector<std::set<int>> graph1(k);
    std::vector<std::set<int>> graph2(k);

    for (int i = 0; i < k; i++)
    {
        if ((static_cast<T>(0.5) * epsilon - vectorDistance(LS[i].begin(), LS[i].end(), LS[0].begin())) < static_cast<T>(0.5) * epsilon)
        {
            graph2[i].insert(i);
            graph1[i].insert(i);

            for (int j = i + 1; j < k; j++)
            {
                // Compute midpoint (m)
                typename std::vector<T>::const_iterator first = LS[i].begin();
                typename std::vector<T>::const_iterator first2 = LS[j].begin();

                for (size_t it = 0; it < m.size(); it++)
                {
                    m[it] = ((*first++) + (*first2++)) / static_cast<T>(2);
                }

                if (vectorDistance(m.begin(), m.end(), LS[0].begin()) < static_cast<T>(0.5) * epsilon)
                {
                    dist = BoundaryProjectedDistance(LS[i].begin(), LS[i].end(), LS[j].begin(), LS[0].begin(), LS[0].end(), static_cast<T>(0.5) * epsilon);
                }
                else
                {
                    dist = static_cast<T>(0.5) * vectorDistance(LS[i].begin(), LS[i].end(), LS[j].begin());
                }

                if (dist < static_cast<T>(0.5) * epsilon)
                {
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

    for (typename std::vector<std::set<int>>::const_iterator it = graph1.begin(); it != graph1.end(); ++it)
    {
        graph3.push_back(std::vector<int>(it->begin(), it->end()));
    }

    for (typename std::vector<std::set<int>>::const_iterator it = graph2.begin(); it != graph2.end(); ++it)
    {
        graph4.push_back(std::vector<int>(it->begin(), it->end()));
    }

    std::list<std::vector<std::vector<int>>> graphs;
    graphs.push_back(graph3);
    graphs.push_back(graph4);
    return graphs;
}
