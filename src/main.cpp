/*
 * hybrid-Lambda is used to simulate gene trees given species network under 
 * coalescent process.
 * 
 * Copyright (C) 2010, 2011, 2012, 2013 Sha (Joe) Zhu
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

#include"param.hpp"
#include"sim_gt.hpp"
#include"freq.hpp"
#include"figure.hpp"
#include"seg-site.hpp"

using namespace hybridLambda;
//using namespace std;



//bool reproduce_GENE_trees;
//string gene_tree_file;

int main(int argc, char *argv[]){
	//remove("seg-sites");

	vector <string> gt_tree_str_s;
	vector <string> mt_tree_str_s;
	
	if (argc==1 ){
		hybridLambda::print_help();
	}	//else, proceed

    try {
	    hybridLambda::param hybrid_para(argc, argv);
	    sim::param sim_para(argc, argv);
	    figure::param figure_para(argc, argv);
	    freq::param freq_para(argc,argv);
	    action_board my_action(argc,argv);
	    
	    if (hybrid_para.print_tree){
			Net new_net_dummy(sim_para.net_str);
			new_net_dummy.print_all_node();
			return EXIT_SUCCESS;
		}
		
		if (hybrid_para.plot_bool){
			figure_para.plot(sim_para.net_str);
			return EXIT_SUCCESS;
		}
		time_t start_time = time(0);
		if (hybrid_para.simulation_bool){
		    sim_n_gt simd_gt_tree_str_s(sim_para,my_action);
			gt_tree_str_s=simd_gt_tree_str_s.gt_string_coal_unit_s;
			mt_tree_str_s=simd_gt_tree_str_s.gt_string_mut_num_s;
		}
		time_t sim_end_time = time(0);
		//std::cout << "Simulation took about " << end_time - start_time  << " second(s)" << std::endl;
				
		if (hybrid_para.freq_bool){
			//frequencies

		}
		time_t freq_end_time=time(0);	
		
		if (hybrid_para.seg_bool){
			//seggreating data were generated
			
		}
		time_t seg_end_time =time(0);
		
		if (hybrid_para.log_bool){          
			std::ofstream log_file;
			log_file.open (hybrid_para.log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
			if (hybrid_para.simulation_bool){
				//simulation 
				//files were saved at where
				log_file << "Simulation took about " << sim_end_time - start_time << " second(s) \n";
			}
			if (hybrid_para.freq_bool){
				//frequencies
				log_file << "Computing topology frequency took about " << freq_end_time - sim_end_time << " second(s) \n";
				
			}
			if (hybrid_para.seg_bool){
				//seggreating data were generated
				log_file << "Generating segregating site data took about " << seg_end_time - freq_end_time << " second(s) \n";
			}
			log_file.close();
			int sys=system("cat log_file");		
		}
    }
    catch (const exception &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}


