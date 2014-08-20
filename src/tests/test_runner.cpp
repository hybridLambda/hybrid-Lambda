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
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>

using namespace CppUnit;

int new_forest_counter = 0;
int delete_forest_counter = 0;


int main(void) {
	TestResult controller;

	TestResultCollector result;
	controller.addListener(&result);

	BriefTestProgressListener progress;
	controller.addListener(&progress);

	TestRunner runner;
	runner.addTest( TestFactoryRegistry::getRegistry().makeTest() );
	runner.run(controller);

	CompilerOutputter outputter(&result, std::cerr);
	outputter.write();

	return result.wasSuccessful() ? 0 : 1;
}

// EOF
