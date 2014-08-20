#include <cppunit/TestCase.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../src/figure.hpp"
#pragma GCC diagnostic ignored "-Wwrite-strings"

class TestFigure : public CppUnit::TestCase {
  
    CPPUNIT_TEST_SUITE( TestFigure );
    
    CPPUNIT_TEST( good_CMD ); 
    CPPUNIT_TEST( bad_CMD ); 

    CPPUNIT_TEST_SUITE_END();

 private:

    void good_CMD() {
        char* argv1[] = { "hybrid-Lambda", "-spcu", "trees/7_tax_sp_nt1_para", "-dotF", "blah", "-branch", "-label" };
        CPPUNIT_ASSERT_NO_THROW ( Figure( 6, argv1 ) ); 
        CPPUNIT_ASSERT_NO_THROW ( Figure( 5, argv1 ) );

        char* argv2[] = { "hybrid-Lambda", "-spcu", "trees/7_tax_sp_nt1_para", "-dotF", "blah", "-label", "-branch" };
        CPPUNIT_ASSERT_NO_THROW ( Figure( 6, argv1 ) ); 
        CPPUNIT_ASSERT_NO_THROW ( Figure( 5, argv1 ) );

        
    }
    
    void bad_CMD(){
        char* argv1[] = { "hybrid-Lambda", "-spcu", "trees/7_tax_sp_nt1_para", "-dotF", "blah", "-branch", "-label" };
        CPPUNIT_ASSERT_THROW ( Figure( 7, argv1 ) , std::invalid_argument ); // too many options
        CPPUNIT_ASSERT_THROW ( Figure( 4, argv1 ) , std::invalid_argument ); // missing file name for -dotF
    }
};

CPPUNIT_TEST_SUITE_REGISTRATION( TestFigure );
