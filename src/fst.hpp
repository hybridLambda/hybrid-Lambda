#include<boost/math/special_functions/binomial.hpp>


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

double ETw( double lambdaA, double lambdaAB, double tau );

double ETb( double lambdaAB, double tau );

double FST_indirect( double lambdaA, double lambdaAB, double tau );

double FST( double lambdaA, double lambdaAB, double tau );

