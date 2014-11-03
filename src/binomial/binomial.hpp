//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SF_BINOMIAL_HPP
#define BOOST_MATH_SF_BINOMIAL_HPP

#ifdef _MSC_VER
#pragma once
#endif

//#include "factorials.hpp"
#include "unchecked_factorial.hpp"
#include "../src/sim_gt.hpp"
//#include "beta.hpp"
//#include <boost/math/policies/error_handling.hpp>

namespace boost{ namespace math{

template <class T>
T binomial_coefficient(T n, T k)
{
   //BOOST_STATIC_ASSERT(!boost::is_integral<T>::value);
   //BOOST_MATH_STD_USING
   static const char* function = "boost::math::binomial_coefficient<%1%>(unsigned, unsigned)";
   if(k > n)
     throw std::invalid_argument ( "boost::math::binomial_coefficient The binomial coefficient is undefined for k > n" );

   T result;
   if((k == 0) || (k == n))
      return 1;
   if((k == 1) || (k == n-1))
      return n;

   if(n <= max_factorial<T>::value)
   {
      // Use fast table lookup:
      result = unchecked_factorial<T>(n);
      result /= unchecked_factorial<T>(n-k);
      result /= unchecked_factorial<T>(k);
   }
   else
   {
      // Use the beta function:
      if(k < n - k)
         //result = k * beta(static_cast<T>(k), static_cast<T>(n-k+1));
         result = k * Beta((double)(k), (double)(n-k+1));
      else
         //result = (n - k) * beta(static_cast<T>(k+1), static_cast<T>(n-k));
         result = (n - k) * Beta((double)(k+1), (double)(n-k));
      if(result == 0)
        throw std::invalid_argument ( "Over flow" );
         //return policies::raise_overflow_error<T>(function, 0, pol);
      result = 1 / result;
   }
   // convert to nearest integer:
   return ceil(result - 0.5f);
}

} // namespace math
} // namespace boost


#endif // BOOST_MATH_SF_BINOMIAL_HPP
