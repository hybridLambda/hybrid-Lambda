/*
 * hybrid-Lambda is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2010 -- 2014 Sha (Joe) Zhu
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

//#include "sim_gt.hpp"
#include"fst.hpp"
/*!
 * Assume two populations A and B have been isolated until time tau in the past as measured from the present. 
 * Assume also that the same coalescent process is operating in populations A and B. 
 * Let TW denote the time until coalescence for two lines when drawn from the same population, 
 * and Tb when drawn from different populations. 
 * Let lambdaA denote the coalescence rate for two lines in population A, and 
 * lambdaAB for the common ancestral population AB. 
 * For the Beta(2 − alpha, alpha)-coalescent, lambdaA = 1, for the point-mass process lambdaA = psi^2. One now obtains
 * ETw exptected value of Tw
 * ETw = (1 - exp(-lambdaA * tau) * lambdaA^{-1} + exp(-lambdaA * tau) * (tau + lambdaAB^{-1})
 * 
 */

double lambda( double alpha ){
    return exp(log(boost::math::binomial_coefficient<double>(unsigned(2),unsigned(2)))+log(Beta(2-alpha,2-2+alpha)) - log(Beta(2.0-alpha,alpha)));
    //return Beta(2-alpha, alpha)/Beta(2.0-alpha,alpha);
}


double ETw( double alphaA, double alphaAB, double tau ){
    //exp(log(boost::math::binomial_coefficient<double>(unsigned(b_i),unsigned(k_i)))+log(Beta(k_i-para,b_i-k_i+para)) - log(Beta(2.0-para,para)))
    double lambdaA = lambda( alphaA );
    double lambdaAB = lambda( alphaAB );
    
    return ( 1 - exp( -lambdaA * tau ) ) / lambdaA + exp( -lambdaA * tau) * ( tau + 1 / lambdaAB);
} 

double ETb( double alphaAB, double tau ){
    double lambdaAB = lambda( alphaAB );
    return tau + 1/ lambdaAB;
}

double FST_indirect( double alphaA, double alphaAB, double tau ){
    double lambdaA = lambda( alphaA );
    double lambdaAB = lambda( alphaAB );
    return ( 1 - ETw( lambdaA, lambdaAB, tau) / ETb( lambdaAB, tau) );
}

double FST( double alphaA, double alphaAB, double tau ){
    double lambdaA = lambda( alphaA );
    double lambdaAB = lambda( alphaAB );
    cout << "lambdaA = " << lambdaA <<endl;
    return ( 1 - exp( -tau ) ) * ( tau / ( 1 + tau ) );
    //return ( 1 - exp(-lambdaA * tau) ) * ( 1 - 1 / ( tau + 1 / lambdaAB ) / lambdaA );
    
}

