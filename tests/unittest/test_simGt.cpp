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

#include "../../src/sim_gt.hpp"

#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestSimGt : public CppUnit::TestCase {
    CPPUNIT_TEST_SUITE( TestSimGt );
    CPPUNIT_TEST ( testInitialization );
    CPPUNIT_TEST ( testLambdaCalculations );
    
    CPPUNIT_TEST_SUITE_END();

  private:
    simTree * testSubjectPtr;
  public:
    void setUp() {
        testSubjectPtr = new simTree();
    }
    
    void tearDown() {
        delete testSubjectPtr;
    }

    void testInitialization(){
    
    }
    
    void testLambdaCalculations(){
        //if p = 1: lambda_{b, b}  =  1,  and all other lambda_{b,k} are 0
        CPPUNIT_ASSERT_EQUAL(1.0, testSubjectPtr->lambdaPsi(2, 2, 1) );
        //CPPUNIT_ASSERT_ASSERTION_FAIL( testSubjectPtr->lambdaPsi(2, 1, 1) ); // This should not happen, k should be always at least 2.
        //CPPUNIT_ASSERT_ASSERTION_FAIL( testSubjectPtr->lambdaPsi(2, 3, 1) );

        CPPUNIT_ASSERT_EQUAL(1.0, testSubjectPtr->lambdaPsi(3, 3, 1) );
        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(3, 2, 1) ); 

        CPPUNIT_ASSERT_EQUAL(1.0, testSubjectPtr->lambdaPsi(30, 30, 1) );
        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(30, 10, 1) );
        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(30, 3, 1) );

        CPPUNIT_ASSERT_EQUAL(1.0, testSubjectPtr->lambdaPsi(100, 100, 1) );
        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(100, 99, 1) );
        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(100, 97, 1) );
        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(100, 93, 1) );

        //if  p = 0:  lambda_{b, 2} =  b*(b-1)/2, and all other lambda_{b,k} are 0
        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(3, 3, 0) );
        CPPUNIT_ASSERT_EQUAL(3, (int)testSubjectPtr->lambdaPsi(3, 2, 0) );

        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(4, 4, 0) );
        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(4, 3, 0) );
        CPPUNIT_ASSERT_EQUAL(6, (int)testSubjectPtr->lambdaPsi(4, 2, 0) );

        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(40, 40, 0) );
        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(40, 39, 0) );
        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(40, 20, 0) );
        CPPUNIT_ASSERT_EQUAL(0.0, testSubjectPtr->lambdaPsi(40, 9, 0) );
        CPPUNIT_ASSERT_EQUAL(780, (int)testSubjectPtr->lambdaPsi(40, 2, 0) );
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestSimGt );


