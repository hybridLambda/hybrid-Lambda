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
#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "src/hybridLambda.hpp"
#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestHybirdLambda : public CppUnit::TestCase {
  
    CPPUNIT_TEST_SUITE( TestHybirdLambda );
    
    CPPUNIT_TEST( test_CMD ); 
    //CPPUNIT_TEST( bad_CMD ); 

    CPPUNIT_TEST_SUITE_END();

 private:

    void test_CMD() {
        char* argv1[] = { "hybrid-Lambda", "-spcu", "../trees/7_tax_sp_nt1_para", "-seg", "-o", "blah", "-sim_mut_unit", "-sim_num_gener" };
        CPPUNIT_ASSERT_NO_THROW ( HybridLambda( 8, argv1 ) ); 
        CPPUNIT_ASSERT_NO_THROW ( HybridLambda( 7, argv1 ) ); 
        CPPUNIT_ASSERT_NO_THROW ( HybridLambda( 6, argv1 ) );
        CPPUNIT_ASSERT_NO_THROW ( HybridLambda( 4, argv1 ) );
        CPPUNIT_ASSERT_NO_THROW ( HybridLambda( 3, argv1 ) );

        
        CPPUNIT_ASSERT_THROW ( HybridLambda( 5, argv1 ) , std::invalid_argument ); 
        CPPUNIT_ASSERT_THROW ( HybridLambda( 2, argv1 ) , std::invalid_argument ); 
        //char* argv2[] = { "hybrid-Lambda", "-spcu", "trees/7_tax_sp_nt1_para", "-dotF", "blah", "-label", "-branch" };
        //CPPUNIT_ASSERT_NO_THROW ( HybridLambda( 6, argv1 ) ); 
        //CPPUNIT_ASSERT_NO_THROW ( HybridLambda( 5, argv1 ) );        
    }
    
    //void bad_CMD(){
        //char* argv1[] = { "hybrid-Lambda", "-spcu", "../trees/7_tax_sp_nt1_para", "-dot", "-o", "blah", "-branch", "-label" };

    //}
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestHybirdLambda );
