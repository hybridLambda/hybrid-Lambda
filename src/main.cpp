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
#include"fst.hpp"

using namespace hybridLambda;

int main(int argc, char *argv[]){

	vector <string> gt_tree_str_s;
	vector <string> mt_tree_str_s;
	
	if ( argc==1 ){
		hybridLambda::print_help();
	}	//else, proceed

    try {
	    hybridLambda::param hybrid_para(argc, argv);
	    figure::param figure_para(argc, argv);
	    freq::param freq_para(argc,argv);
	    action_board my_action(argc,argv);
	    seg::param seg_para(argc,argv);
        double Fst;

   	    if (hybrid_para.seg_bool){my_action.sim_num_mut_bool=true;}

		time_t start_time = time(0);
		if (hybrid_para.simulation_bool){
			sim::param sim_para(argc, argv);	

			Net new_net_dummy(sim_para.net_str);
			if (hybrid_para.print_tree){
				new_net_dummy.print_all_node();
				return EXIT_SUCCESS;
			}
			
			if (hybrid_para.plot_bool){
				figure_para.plot(sim_para.net_str);
				return EXIT_SUCCESS;
			}
			
			if (!new_net_dummy.is_ultrametric){
				cout<<"WARNING! NOT ULTRAMETRIC!!!"<<endl;
				//return EXIT_SUCCESS;
			}
			
		    sim_n_gt simd_gt_tree_str_s(sim_para,my_action);
			gt_tree_str_s = simd_gt_tree_str_s.gt_string_coal_unit_s;
			mt_tree_str_s = simd_gt_tree_str_s.gt_string_mut_num_s;
			
			
			if (my_action.mono_bool  &&  sim_para.sample_size.size()==2){
				cout<<"   A mono     B mono Recip mono     A para     B para  Polyphyly"<<endl;
				for (unsigned int mono_i=0;mono_i<simd_gt_tree_str_s.monophyly.size();mono_i++){
					cout<<setw(9)<<simd_gt_tree_str_s.monophyly[mono_i]<<"  ";
				}
				cout<<endl;
			}
		}
		time_t sim_end_time = time(0);
		
		if (hybrid_para.read_GENE_trees){
			gt_tree_str_s=read_input_lines(hybrid_para.gt_file_name.c_str());
		}
		
		if (hybrid_para.tmrca_bool){
			remove(hybrid_para.tmrca_NAME.c_str());    
			std::ofstream tmrca_file;
			tmrca_file.open (hybrid_para.tmrca_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
			for (size_t i=0;i<gt_tree_str_s.size();i++){
				Net gt(gt_tree_str_s[i]);
				tmrca_file<< gt.Net_nodes.back().absolute_time<<endl;
			}
			tmrca_file.close();
		}
		
		if (hybrid_para.bl_bool){
			remove(hybrid_para.bl_NAME.c_str());    
			std::ofstream bl_file;
			bl_file.open (hybrid_para.bl_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
			for (size_t i=0;i<gt_tree_str_s.size();i++){
				Net gt(gt_tree_str_s[i]);
				double totalbl=0;
				for (size_t node_i = 0 ; node_i < gt.Net_nodes.size(); node_i++){
					totalbl = totalbl + gt.Net_nodes[node_i].brchlen1;
				}
				bl_file << totalbl << endl;
			}
			bl_file.close();
		}
		
		if (hybrid_para.read_mt_trees){
			mt_tree_str_s=read_input_lines(hybrid_para.mt_file_name.c_str());
		}
				
		if (hybrid_para.freq_bool){
			//frequencies
			compute_gt_frequencies(gt_tree_str_s, freq_para.freq_file_name);

		}
		time_t freq_end_time=time(0);	
		if (hybrid_para.seg_bool){
			//seggreating data were generated
			seg_para.create_site_data_dir(mt_tree_str_s);
		}
		time_t seg_end_time =time(0);
        
        if (hybrid_para.fst_bool){
            sim::param sim_para(argc, argv);	
            Net coal_unit_net(sim_para.net_str);
            double tau = coal_unit_net.Net_nodes[0].brchlen1;
            
            //Net para_net(sim_para.para_string);
            //double lambdaA  = para_net.Net_nodes[0].brchlen1;
            //double lambdaAB = para_net.Net_nodes.back().brchlen1;
            //cout<<"tau = " << tau <<endl;
            //cout << "( 1 - exp( -tau ) )= "<<( 1 - exp( -tau ) )<<endl;
            Fst = ( 1 - exp( -tau ) ) * ( tau / ( 1 + tau ) );
            cout << "Expected[Fst] = " << Fst << endl;
        }
		
		if (hybrid_para.log_bool){      
			remove(hybrid_para.log_NAME.c_str());    
			std::ofstream log_file;
			log_file.open (hybrid_para.log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
			if (hybrid_para.simulation_bool){
				//simulation 
				//files were saved at where
				log_file << "Simulation took about " << sim_end_time - start_time << " second(s) \n";
				log_file << "Random seed is :"<<hybrid_para.seed << "\n";
				if ( !hybrid_para.mm_bool ){
					log_file<<"Default Kingman coalescent on all branches. \n";
					}
				if ( !hybrid_para.pop_bool ){
					log_file << "Default population size of 10000 on all branches. \n";
					}
				log_file << "Produced gene tree files: \n";	
				log_file << my_action.gene_tree_file<<"_coal_unit\n";
				
				if (my_action.sim_mut_unit_bool){
					log_file << my_action.gene_tree_file<<"_mut_unit\n";
				}
				if (my_action.sim_num_gener_bool){
					log_file << my_action.gene_tree_file<<"_numb_gener\n";
				}
				if (my_action.sim_num_mut_bool){
					log_file << my_action.gene_tree_file<<"_num_mut\n";
				}
				
			}
			
			if (hybrid_para.tmrca_bool){
				log_file << "TMRCA file is saved at: "<<hybrid_para.tmrca_NAME<<"\n";
			}

			if (hybrid_para.bl_bool){
				log_file << "Total branch length file is saved at: "<<hybrid_para.bl_NAME<<"\n";
			}

			
			if (hybrid_para.freq_bool){
				//frequencies
				log_file << "Computing topology frequency took about " << freq_end_time - sim_end_time << " second(s) \n";
				log_file << "Frequency file is saved at: "<<freq_para.freq_file_name<<"\n";
			}
			if (hybrid_para.seg_bool){
				//seggreating data were generated
				log_file << "Generating segregating site data took about " << seg_end_time - freq_end_time << " second(s) \n";
				log_file << "Segregating site data saved at: "<<seg_para.seg_dir_name<<"\n";
			}
            if (hybrid_para.fst_bool){
                log_file << "Fst = " << Fst << "\n";
            }
            
			log_file.close();
			string showlog="cat "+ hybrid_para.log_NAME;
			int sys=system(showlog.c_str());
		}
    }
    catch (const exception &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}


