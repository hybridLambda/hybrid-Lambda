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

/*! \file main.cpp
 *  \brief Main function of hybrid-Lambda */

#include "hybridLambda.hpp"
#include "fst.hpp"

int main(int argc, char *argv[]){

	if ( argc == 1 ) print_help(); 	//else, proceed

    try {
	    HybridLambda run_hybridLambda ( argc, argv );
        double Fst;

        run_hybridLambda.HybridLambda_core( );
        run_hybridLambda.create_site_data_dir(); // segregating site data were generated				

        run_hybridLambda.extract_tmrca ();
        run_hybridLambda.extract_bl ();
        run_hybridLambda.extract_firstcoal();
        run_hybridLambda.extract_frequency();


/// need to work         
//if (run_hybridLambda.fst_bool){
////sim::param sim_para(argc, argv);	
//Net coal_unit_net(sim_para.net_str);
//double tau = coal_unit_net.NodeContainer[0].brchlen1;

////Net para_net(sim_para.para_string);
////double lambdaA  = para_net.NodeContainer[0].brchlen1;
////double lambdaAB = para_net.NodeContainer.back().brchlen1;
////cout<<"tau = " << tau <<endl;
////cout << "( 1 - exp( -tau ) )= "<<( 1 - exp( -tau ) )<<endl;
//Fst = ( 1 - exp( -tau ) ) * ( tau / ( 1 + tau ) );
//cout << "Expected[Fst] = " << Fst << endl;
//}

		//}			
    }
    catch (const exception &e) {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}

