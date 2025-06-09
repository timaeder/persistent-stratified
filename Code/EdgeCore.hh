#pragma once
#include <vector>
#include <algorithm>
#include <set>
#include "RipsComplex.hh"
#include "unordered_map.hpp"

using namespace std;

bool EdgeDomination
( 
    vector<int>& edge,
    ska::unordered_map< int, vector<int> >& nbh
)
{

    int v0 = edge[0];
    int v1 = edge[1];

    vector<int> Nedge;
    
    std::set_intersection
    (
        nbh[v0].begin(), nbh[v0].end(),
        nbh[v1].begin(), nbh[v1].end(),
        back_inserter(Nedge)
    );


    for (const auto& N : nbh)
    {
        if (N.first != v0 && N.first != v1)
        {
            if (includes((N.second).begin(),(N.second).end(),Nedge.begin(),Nedge.end()))
            {
                return true;
            }

        }
    }

    return false;
}


list<vector< vector<int>>> CplxNbh
(
    const vector< vector<int> >& E,
    const set<int>& V 
)
{
    list<vector< vector<int>>> graph;
    vector< vector<int> > graph1(V.size());
    vector< vector<int> > graph2(V.size());

    for (const int& v : V)
    {
        //graph1[v] = vector<int>{v};
        graph2[v] = vector<int>{v};
    }

    for (const auto& e : E)
    {
        graph2[e[0]].push_back(e[1]);
        graph2[e[1]].push_back(e[0]);
    }

    for (int i = 0; i < graph2.size(); i++)
    {
        sort(graph2[i].begin(), graph2[i].end());
        vector<int>::iterator find_it =
        find(graph2[i].begin(), graph2[i].end(),i);
        copy(find_it, graph2[i].end(),
         back_inserter(graph1[i]));
    }

    graph.push_back(graph1);
    graph.push_back(graph2);


    return graph;
}

list<vector< vector<int>>> RelCplxNbh
(
    const vector< vector<int> >& E,
    const set<int>& V 
)
{
    list<vector< vector<int>>> graph;
    int m = *V.rbegin();
    vector< vector<int> > graph1(m+1);
    vector< vector<int> > graph2(m+1);

    for (const int& v : V)
    {
        //graph1[v] = vector<int>{v};
        graph2[v] = vector<int>{v};
    }

    for (const auto& e : E)
    {
        graph2[e[0]].push_back(e[1]);
        graph2[e[1]].push_back(e[0]);
    }

    for (int i = 0; i < graph2.size(); i++)
    {
        if(graph2[i].size() > 1)
        {
            sort(graph2[i].begin(), graph2[i].end());
            vector<int>::iterator find_it =
            find(graph2[i].begin(), graph2[i].end(),i);
            copy(find_it, graph2[i].end(),
            back_inserter(graph1[i]));
        }
    }

    graph.push_back(graph1);
    graph.push_back(graph2);


    return graph;
}


ska::unordered_map< int, vector<int> > GraphNbh
(
    const vector< vector<int> >& E,
    const set<int>& V
)
{
    ska::unordered_map< int, vector<int> > graph(V.size());

    for (const int& v : V)
    {
        graph[v] = vector<int>{v};
    }

    for (const auto& e : E)
    {
        graph[e[0]].push_back(e[1]);
        graph[e[1]].push_back(e[0]);
    }

    for (auto& N : graph)
    {
        sort(N.second.begin(),N.second.end());
    }


    return graph;
}


ska::unordered_set< vector<int>, VectorHasher > EdgeNbh
(
    vector< int >& edge, 
    ska::unordered_map< int, vector<int> >& nbh
)
{
    int v0 = edge[0];
    int v1 = edge[1];

    ska::unordered_set< vector<int>, VectorHasher > EN;
    vector<int> Nedge;

    set_intersection
    (
        nbh[v0].begin(), nbh[v0].end(),
        nbh[v1].begin(), nbh[v1].end(),
        back_inserter(Nedge)
    );    

    for (const auto& y : Nedge)
    {
        EN.insert({v0,y});
        EN.insert({v1,y});
        EN.insert({y,v0});
        EN.insert({y,v1});
    }

    return EN;

}

set<int> VertexSet
(
    vector< vector<int> >& E
)
{
    set<int> V;

    for (const auto& e : E)
    {
        V.insert(e[0]);
        V.insert(e[1]);
    }
    return V;
}

ska::unordered_map< vector<int>, int, VectorHasher > EdgeFlagCore
(
    vector< vector<int> >& E
)
{
    ska::unordered_map< vector<int>, int, VectorHasher > EC;
    set<int> V = VertexSet(E);

    ska::unordered_map< int, vector<int> > N(V.size());

    for (const int& v : V)
    {
        N[v] = vector<int>{v};
    }

    vector <int> ei(2);

    int n = E.size();

    vector< vector<int> > Ei;

    for (int i = 0; i < n; i++)
    {
        ei = E[i];
        Ei.push_back(ei);

        N[ei[0]].push_back(ei[1]);
        sort(N[ei[0]].begin(),N[ei[0]].end());
        N[ei[1]].push_back(ei[0]);
        sort(N[ei[1]].begin(),N[ei[1]].end());

        ska::unordered_map< int, vector<int> > Ni = N;

        if(!(EdgeDomination(ei,Ni)))
        {
            EC[ei] = i;

            ska::unordered_set< vector<int>, VectorHasher > EN = EdgeNbh(ei,Ni);

            vector< vector<int> >::iterator simp_it = Ei.end()-2;

            while(simp_it != Ei.begin()-1)
            {

                if(EC.find(*simp_it) == EC.end())
                {
                    if(EN.find(*simp_it) != EN.end())
                    {
                        if(!(EdgeDomination(*simp_it,Ni)))
                        {
                            EC[*simp_it] = i;

                            vector<int> Nedge;

                            set_intersection
                            (
                            Ni[simp_it->front()].begin(), Ni[simp_it->front()].end(),
                            Ni[simp_it->back()].begin(), Ni[simp_it->back()].end(),
                            back_inserter(Nedge)
                            );    

                            for (const auto& y : Nedge)
                            {
                                EN.insert({simp_it->front(),y});
                                EN.insert({simp_it->back(),y});
                                EN.insert({y,simp_it->front()});
                                EN.insert({y,simp_it->back()});
                            }
                        }
                        else
                        {
                            vector<int>::iterator find_it = 
                            find(Ni[simp_it->front()].begin(),Ni[simp_it->front()].end(),simp_it->back());
                            Ni[simp_it->front()].erase(find_it);
                            find_it = 
                            find(Ni[simp_it->back()].begin(),Ni[simp_it->back()].end(),simp_it->front());
                            Ni[simp_it->back()].erase(find_it);
                        }
                    }
                    else
                    {
                        vector<int>::iterator find_it = 
                        find(Ni[simp_it->front()].begin(),Ni[simp_it->front()].end(),simp_it->back());
                        Ni[simp_it->front()].erase(find_it);
                        find_it = 
                        find(Ni[simp_it->back()].begin(),Ni[simp_it->back()].end(),simp_it->front());
                        Ni[simp_it->back()].erase(find_it);
                    }
                }
                simp_it--;
            }

        }

    }

    return EC;
}
