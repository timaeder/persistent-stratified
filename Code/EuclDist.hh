#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <type_traits> // For std::remove_const


template<class Iter_T1, class Iter_T2>
typename std::remove_const<typename std::iterator_traits<Iter_T1>::value_type>::type vectorDistance(
    Iter_T1 first,
    Iter_T1 last,
    Iter_T2 first2
) {
    using T = typename std::remove_const<typename std::iterator_traits<Iter_T1>::value_type>::type;

    T ret = T{0};
    while (first != last) {
        T dist = ((*first++) - (*first2++));
        ret += (dist * dist);
    }
    return ret > T{0} ? sqrt(ret) : T{0};
}

template <typename Iter_T>
typename std::remove_const<typename std::iterator_traits<Iter_T>::value_type>::type BoundaryProjectedDistance(
    Iter_T x_begin, Iter_T x_end,
    Iter_T y_begin,
    Iter_T z_begin, Iter_T z_end,
    const typename std::iterator_traits<Iter_T>::value_type& r
) {
    using T = typename std::remove_const<typename std::iterator_traits<Iter_T>::value_type>::type;

    T d = T{0};

    // Compute dot product of (z - m) and (x - y) on the fly
    T ret = T{0};
    T dxy = T{0};
    T dmz = T{0};
    T doz = T{0};

    Iter_T x_it = x_begin;
    Iter_T y_it = y_begin;
    Iter_T z_it = z_begin;

    while (x_it != x_end && z_it != z_end) {
        T m = (*x_it + *y_it) / 2; // Midpoint
        T xmy = *x_it - *y_it;     // Difference vector (x - y)
        T zmm = *z_it - m;         // Difference vector (z - m)

        ret += zmm * xmy;          // Dot product of (z - m) and (x - y)
        dxy += xmy * xmy;          // Distance squared between x and y
        dmz += zmm * zmm;          // Distance squared between z and m

        ++x_it;
        ++y_it;
        ++z_it;
    }

    dxy = dxy > T{0} ? sqrt(dxy) : T{0}; // Final distance between x and y
    doz = (dxy != T{0}) ? std::fabs(ret / dxy) : T{0}; // Projection distance
    dmz = dmz > T{0} ? sqrt(dmz) : T{0}; // Final distance between z and m

    // Ensure non-negative terms for square roots
    T term1 = (r * r - doz * doz > T{0}) ? sqrt(r * r - doz * doz) : T{0};
    T term2 = (dmz * dmz - doz * doz > T{0}) ? sqrt(dmz * dmz - doz * doz) : T{0};

    // Final computation
    d = (term1 - term2) * (term1 - term2) + dxy * dxy / 4;

    // Return result
    return d > T{0} ? sqrt(d) : T{0};
}
