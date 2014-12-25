#ifndef hybridLambda_src_macros
#define hybridLambda_src_macros

/* Used when building an R package (RBUILD) */
#ifdef RBUILD

// Include Rcpp Headers for Rcout.
#include "Rcpp.h"

// Suppress debug output.
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && Rcpp::Rcout

// Assure that assertions are deactivated.
#ifndef NDEBUG
#define NDEBUG
#endif

#else
/* Used for normal compilation for hybridLambda */

// Unless compiled with options NDEBUG, we will produce a debug output using 
// 'dout' instead of cout and execute (expensive) assert statements.
#ifndef NDEBUG

// Debug mode
#ifdef UNITTEST // No debug output in unittests
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#else           // Produce debug output
#define dout std::cout
#endif

#else
// Normal Mode
#pragma GCC diagnostic ignored "-Wunused-value"
#define dout 0 && std::cout
#endif

#endif

#endif
