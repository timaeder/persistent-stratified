#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

template<class Iter_T, class Iter2_T>
long double vectorDistance(Iter_T first, Iter_T last, Iter2_T first2) {
  long double ret = 0.0;
  while (first != last) {
    long double dist = ((*first++) - (*first2++));
    ret += (dist * dist);
  }
  return ret > 0.0 ? sqrtl(ret) : 0.0;
}

template<class Iter_T, class Iter2_T>
long double vectorOneDistance(Iter_T first, Iter_T last, Iter2_T first2) {
  long double ret = 0.0;
  while (first != last) {
    long double dist = (*first++) - (*first2++);
    ret += std::fabs(dist);
  }
  return ret > 0.0 ? ret : 0.0;
}

template<class Iter_T, class Iter2_T>
long double vectorSupDistance(Iter_T first, Iter_T last, Iter2_T first2) {
  long double ret = 0.0;
  while (first != last) {
    long double dist = std::fabs((*first++) - (*first2++));
    if (ret < dist)
    {
      ret = dist;
    }
  }
  return ret > 0.0 ? ret : 0.0;
}

long double MinRad(std::vector<std::vector<long double>>& S)
{
  std::vector<long double> dist;
  std::vector<long double> max;
  int n = S.size();
  for( int i = 0; i < n; i++)
  {
    for(int j = 0; j < n; j++)
    {
      dist.push_back(vectorDistance(S[i].begin(),S[i].end(),S[j].begin()));
    }
    max.push_back(*std::max_element(dist.begin(),dist.end()));
    dist.clear();
  }
return *std::min_element(max.begin(),max.end());
}

template<class T, size_t N>
constexpr size_t size(T (&)[N]) { return N; }


long double MaxDist(std::vector< std::vector<long  double> >& S, std::vector<long double> p)
{
  long double dist = 0.0;
  long double max = 0.0;
  for (int i = 0; i < S.size(); i++)
  {
    dist = vectorDistance(p.begin(),p.end(),S[i].begin());
    if (dist > max)
    {
      max = dist;
    }
    
  }
  return max;
}

template<class T> 
class midVec
{
public:
    T operator()(T &lhs, const T &rhs) 
    {
    return (lhs + rhs)/2;
    }
};

template <typename Iter_T>
auto G(
    Iter_T x_begin, Iter_T x_end,
    Iter_T y_begin,
    Iter_T z_begin, Iter_T z_end,
    const typename std::iterator_traits<Iter_T>::value_type& r
) -> typename std::iterator_traits<Iter_T>::value_type
{
    using T = typename std::iterator_traits<Iter_T>::value_type;

    T d = 0.0;
    std::vector<T> m;
    std::vector<T> xmy;

    // Compute midpoint (m) and difference vector (xmy)
    Iter_T x_it = x_begin;
    Iter_T y_it = y_begin;
    while (x_it != x_end)
    {
        m.push_back((*x_it + *y_it) / 2);
        xmy.push_back(*x_it - *y_it);
        ++x_it;
        ++y_it;
    }

    // Compute z - m (zmm)
    std::vector<T> zmm;
    Iter_T z_it = z_begin;
    auto m_it = m.begin();
    while (z_it != z_end)
    {
        zmm.push_back(*z_it - *m_it);
        ++z_it;
        ++m_it;
    }

    // Compute dot product of zmm and xmy
    T ret = 0.0;
    for (size_t i = 0; i < zmm.size(); ++i)
    {
        ret += zmm[i] * xmy[i];
    }

    // Compute distances
    T dxy = vectorDistance(x_begin, x_end, y_begin);
    T doz = (dxy != 0.0) ? std::fabs(ret / dxy) : 0.0;
    T dmz = vectorDistance(m.begin(), m.end(), z_begin);

    // Ensure non-negative terms for square roots
    T term1 = (r * r - doz * doz > 0.0) ? sqrt(r * r - doz * doz) : 0.0;
    T term2 = (dmz * dmz - doz * doz > 0.0) ? sqrt(dmz * dmz - doz * doz) : 0.0;

    // Final computation
    d = (term1 - term2) * (term1 - term2) + dxy * dxy / 4;

    // Return result
    return d > 0.0 ? sqrt(d) : 0.0;
}

template <typename Iter_T>
auto BoundaryProjectedDistance(
    Iter_T x_begin, Iter_T x_end,
    Iter_T y_begin,
    Iter_T z_begin, Iter_T z_end,
    const typename std::iterator_traits<Iter_T>::value_type& r
) -> typename std::iterator_traits<Iter_T>::value_type
{
    using T = typename std::iterator_traits<Iter_T>::value_type;

    T d = 0.0;
    std::vector<T> m;
    std::vector<T> xmy;

    // Compute midpoint (m) and difference vector (xmy)
    Iter_T x_it = x_begin;
    Iter_T y_it = y_begin;
    while (x_it != x_end)
    {
        m.push_back((*x_it + *y_it) / 2);
        xmy.push_back(*x_it - *y_it);
        ++x_it;
        ++y_it;
    }

    // Compute z - m (zmm)
    std::vector<T> zmm;
    Iter_T z_it = z_begin;
    auto m_it = m.begin();
    while (z_it != z_end)
    {
        zmm.push_back(*z_it - *m_it);
        ++z_it;
        ++m_it;
    }

    // Compute dot product of zmm and xmy
    T ret = 0.0;
    for (size_t i = 0; i < zmm.size(); ++i)
    {
        ret += zmm[i] * xmy[i];
    }

    // Compute distances
    T dxy = vectorDistance(x_begin, x_end, y_begin);
    T doz = (dxy != 0.0) ? std::fabs(ret / dxy) : 0.0;
    T dmz = vectorDistance(m.begin(), m.end(), z_begin);

    // Ensure non-negative terms for square roots
    T term1 = (r * r - doz * doz > 0.0) ? sqrt(r * r - doz * doz) : 0.0;
    T term2 = (dmz * dmz - doz * doz > 0.0) ? sqrt(dmz * dmz - doz * doz) : 0.0;

    // Final computation
    d = (term1 - term2) * (term1 - term2) + dxy * dxy / 4;

    // Return result
    return d > 0.0 ? sqrt(d) : 0.0;
}

long double Hausdorff
(
  std::vector< std::vector<long double> >& S,
  std::vector< std::vector<long double> >& T
)
{

  long double dist = 0.0;
  long double max = 0.0;
  long double min = 1000.0;
  long double M;

  for (int i = 0; i < S.size(); i++)
  {
    for (int j = 0; j < T.size(); j++)
    {
      dist = vectorDistance(S[i].begin(),S[i].end(),T[j].begin());
      if(dist < min)
      {
        min = dist;
      }
    }
    if(max < min)
    {
      max = min;
    }
  }

  M = max;
  max = 0.0;
  min = 1000.0;
  
  for (int i = 0; i < T.size(); i++)
  {
    for (int j = 0; j < S.size(); j++)
    {
      dist = vectorDistance(T[i].begin(),T[i].end(),S[j].begin());
      if(dist < min)
      {
        min = dist;
      }
    }
    if(max < min)
    {
      max = min;
    }
  }


  if (M < max)
  {
    return max;
  }
  else
  {
    return M;
  }

}