#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <numeric>

/*
vectorDistance is more or less taken from:

C++ Cookbook
by D. Ryan Stephens, Christopher Diggins, Jonathan Turkanis, Jeff Cogswell
Released November 2005
Publisher(s): O'Reilly Media, Inc.
ISBN: 9780596007614

*/


template<class Iter> typename std::iterator_traits<Iter>::value_type 
vectorDistance(Iter first, Iter last, Iter first2) {
  typename std::iterator_traits<Iter>::value_type ret = 0.0;
  while (first != last) {
    typename std::iterator_traits<Iter>::value_type dist = ((*first++) - (*first2++));
    ret += (dist * dist);
  }
  return ret > 0.0 ? sqrt(ret) : 0.0;
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


//Function to correctly determine the simplex diameter as in "Approximating Local Homology from Samples" by Skraba and Wang.
template<class Iter1> typename std::iterator_traits<Iter1>::value_type 
G
(
  Iter1 first1,
  Iter1 last1,
  Iter1 first2,
  Iter1 first3,
  Iter1 last3,
  const typename std::iterator_traits<Iter1>::value_type r
)
{
  std::vector<typename std::iterator_traits<Iter1>::value_type> m;
  std::vector<typename std::iterator_traits<Iter1>::value_type> xmy;
  std::vector<typename std::iterator_traits<Iter1>::value_type> zmm;

  std::transform(first1, last1, first2, std::back_inserter(m), midVec<const typename std::iterator_traits<Iter1>::value_type>());

  std::transform(first1, last1, first2, std::back_inserter(xmy), std::minus<const typename std::iterator_traits<Iter1>::value_type>());

  std::transform(first3, last3, m.begin(), std::back_inserter(zmm), std::minus<const typename std::iterator_traits<Iter1>::value_type>());

  Iter1 first = zmm.begin();
  Iter1 last = zmm.end();
  Iter1 first_2 = xmy.begin();

  const typename std::iterator_traits<Iter1>::value_type ret
   = std::inner_product(zmm.begin(), zmm.end(), xmy.begin(), 0);
  
  const typename std::iterator_traits<Iter1>::value_type dxy = vectorDistance(first1,last1,first2);

  const typename std::iterator_traits<Iter1>::value_type doz = std::abs(ret/dxy);
  
  Iter1 Iter1_m = m.begin();
  Iter1 Iter2_m = m.end();

  const typename std::iterator_traits<Iter1>::value_type dmz = vectorDistance(Iter1_m,Iter2_m,first3);
  
  const typename std::iterator_traits<Iter1>::value_type d
   = (sqrt(r*r - doz*doz) - sqrt(dmz*dmz - doz*doz))
      *(sqrt(r*r - doz*doz) - sqrt(dmz*dmz - doz*doz)) 
      + dxy*dxy/4;

  return d > 0.0 ? sqrt(d) : 0.0;
}
