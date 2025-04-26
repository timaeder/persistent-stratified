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

    //in order to "reduce homology" ?
    LocCopairs[0].erase(LocCopairs[0].end()-1);

    return LocCopairs;
}
