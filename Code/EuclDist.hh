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
