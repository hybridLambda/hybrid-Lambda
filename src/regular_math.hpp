/*
 * hybrid-Lambda is used to simulate gene trees given species network under
 * coalescent process.
 *
 * Copyright (C) 2010 -- 2015 Sha (Joe) Zhu
 *
 * This file is part of hybrid-Lambda.
 *
 * hybrid-Lambda is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#include "global.hpp"

/*! \brief Compute factorial of a \return double a! */
template < class T > T factorial ( T a ){
    if (a > 1) return (a * factorial (a-1));
    else       return (1);
}

/*! \brief Compute a permutations of n \return double */
template < class T > T n_permu_a ( T n, T a ){
    if   ( a > 1 ) return (n*n_permu_a(n-1,a-1));
    else if (a==1) return (n);
    else           return (1);
}

/*! \brief Compute n choose k \return double */
template < class T > T n_choose_k ( T n, T k ){
    if ( k < ( n/2 ) ) return (n_choose_k(n,n-k));
    else               return (n_permu_a(n,k)/factorial(k));
}
