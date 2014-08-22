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

#include"hybridLambda.hpp"
#include"freq.hpp"
#include"figure.hpp"
#include"seg-site.hpp"
#include"fst.hpp"


int main(int argc, char *argv[]){

	if ( argc == 1 ) print_help(); 	//else, proceed

    try {
	    HybridLambda hybrid_para ( argc, argv );
        //hybrid_para.parse();
        Figure figure_para ( argc, argv );
	    Freq freq_para ( argc, argv );
	    seg::param seg_para(argc,argv);
        double Fst;

   	    //if ( hybrid_para.seg_bool ){my_action.sim_num_mut_bool=true;}

		//time_t start_time = time(0);
		if (hybrid_para.simulation_bool){
			//sim::param sim_para(argc, argv);	

			Net new_net_dummy(hybrid_para.parameters()->net_str);
			if (hybrid_para.print_tree){
				new_net_dummy.print_all_node();
				return EXIT_SUCCESS;
			}
			
			if (hybrid_para.plot_bool){
				figure_para.plot(hybrid_para.parameters()->net_str);
				return EXIT_SUCCESS;
			}
			
			if (!new_net_dummy.is_ultrametric){
				cout<<"WARNING! NOT ULTRAMETRIC!!!"<<endl;
				//return EXIT_SUCCESS;
			}
			
            hybrid_para.HybridLambda_core( );
			
			if (hybrid_para.simulation_jobs()->mono() ){
                if ( hybrid_para.parameters()->sample_size.size() == 2 ){ //\todo check population structure is a species tree
                    cout << "   A mono     B mono Recip mono     A para     B para  Polyphyly" << endl;
                    for ( size_t mono_i = 0; mono_i < hybrid_para.monophyly.size(); mono_i++){
                        cout << setw(9) << hybrid_para.monophyly[mono_i] << "  ";
                    }
                    cout << endl;
                }
                else {
                    throw std::invalid_argument(std::string(" -mono flag can only apply to species tree of two population") );
                    }
			}
		}
		//time_t sim_end_time = time(0);
		
		if (hybrid_para.read_GENE_trees) hybrid_para.gt_tree_str_s = read_input_lines(hybrid_para.gt_file_name.c_str());
		
        hybrid_para.extract_tmrca ();
        hybrid_para.extract_bl ();
        hybrid_para.extract_firstcoal();
        
        //frequencies			
		if ( hybrid_para.freq_bool ) freq_para.compute_gt_frequencies( hybrid_para.gt_tree_str_s );       

		if (hybrid_para.read_mt_trees) hybrid_para.mt_tree_str_s = read_input_lines(hybrid_para.mt_file_name.c_str());
        
		//time_t freq_end_time=time(0);	
		if (hybrid_para.seg_bool){ 	//seggreating data were generated
			//seg_para.create_site_data_dir(mt_tree_str_s);
            seg_para.create_site_data_dir( hybrid_para.mt_tree_str_s );
		}
		//time_t seg_end_time =time(0);


/// need to work         
//if (hybrid_para.fst_bool){
////sim::param sim_para(argc, argv);	
//Net coal_unit_net(sim_para.net_str);
//double tau = coal_unit_net.Net_nodes[0].brchlen1;

////Net para_net(sim_para.para_string);
////double lambdaA  = para_net.Net_nodes[0].brchlen1;
////double lambdaAB = para_net.Net_nodes.back().brchlen1;
////cout<<"tau = " << tau <<endl;
////cout << "( 1 - exp( -tau ) )= "<<( 1 - exp( -tau ) )<<endl;
//Fst = ( 1 - exp( -tau ) ) * ( tau / ( 1 + tau ) );
//cout << "Expected[Fst] = " << Fst << endl;
//}

        //if (hybrid_para.log_bool){      
			//remove(hybrid_para.log_NAME.c_str());    
			//std::ofstream log_file;
			//log_file.open (hybrid_para.log_NAME.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
			//if (hybrid_para.simulation_bool){
				////simulation 
				////files were saved at where
				//log_file << "Simulation took about " << sim_end_time - start_time << " second(s) \n";
				//log_file << "Random seed is :"<<hybrid_para.seed << "\n";
				//if ( !hybrid_para.mm_bool ){
					//log_file<<"Default Kingman coalescent on all branches. \n";
					//}
				//if ( !hybrid_para.pop_bool ){
					//log_file << "Default population size of 10000 on all branches. \n";
					//}
				
			//}
						
			//if (hybrid_para.seg_bool){
				////seggreating data were generated
				//log_file << "Generating segregating site data took about " << seg_end_time - freq_end_time << " second(s) \n";
				//log_file << "Segregating site data saved at: "<<seg_para.seg_dir_name<<"\n";
			//}
            //if (hybrid_para.fst_bool){
                //log_file << "Fst = " << Fst << "\n";
            //}
            
			//log_file.close();
			//string showlog="cat "+ hybrid_para.log_NAME;
			//int sys=system(showlog.c_str());
		//}
    }
    catch (const exception &e)
    {
      std::cerr << "Error: " << e.what() << std::endl;
      return EXIT_FAILURE;
    }
}

void print_example(){
	cout<<"Examples:"<<endl;
	cout<<""<<endl;	
	cout<<"hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 3 -seed 2 -o example1"<<endl;	
	cout<<"hybrid-Lambda -spcu trees/4_tax_sp_nt1_para -o example2 -num 2 -mu 0.00003 -sim_mut_unit -sim_num_mut"<<endl;
	cout<<"hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 100 -pop 25000 -sim_num_gener"<<endl;
	cout<<"hybrid-Lambda -spng '(A:50000,B:50000)r;' -pop '(A:50000,B:50000)r:40000;'"<<endl;
	cout<<"hybrid-Lambda -spcu '((((A:1.1,B:1.1):2.1,a:2.2):1.1,13D:.2):.3,4:.3);' -S 2 4 3 6 5"<<endl;
	cout<<"hybrid-Lambda -spcu '(A:1,B:1)r;' -mm '(A:1.9,B:.2)r:2;' -S 3 4"<<endl;
	cout<<"hybrid-Lambda -spcu trees/7_tax_sp_nt1_para -dot -branch"<<endl;	
	cout<<"hybrid-Lambda -spcu trees/4_tax_sp1.tre -num 1000 -o GENE_TREE_FILE -f"<<endl;	
	cout<<"hybrid-Lambda -spcu trees/4_tax_sp1.tre -num 1000 -o GENE_TREE_FILE -fF FRENQUENCY_FILE"<<endl;	
	cout<<"hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 1000 -o GENE -fF OUTPUT"<<endl;	
	cout<<"hybrid-Lambda -gt GENE_coal_unit -f "<<endl;	
	cout<<"hybrid-Lambda -spcu '(A:5,B:5)r;' -mono -num 100 -mm .1 -S 4 4"<<endl;
	cout<<endl;
}

void print_option(){
	cout<<setw(20)<<"-h or -help"<<"  --  "<<"Help. List the following content."<<endl;
	cout<<setw(20)<<"-spcu INPUT"<<"  --  "<<"Input the species network/tree string through command line or a file."<<endl;
	cout<<setw(26)<<" "<<"Branch lengths of the INPUT are in coalescent unit."<<endl;
	cout<<setw(20)<<"-spng INPUT"<<"  --  "<<"Input the species network/tree string through command line or a file. "<<endl;
	cout<<setw(26)<<" "<<"Branch lengths of the INPUT are in number of generation."<<endl;
	cout<<setw(20)<<"-pop INPUT"<<"  --  "<<"Population sizes are defined by a single numerical constant, "<<endl;
	cout<<setw(26)<<" "<<"or a string which specifies the population size on each branch. "<<endl;
	cout<<setw(26)<<" "<<"The string can be input through command line or a file. "<<endl;
	cout<<setw(26)<<" "<<"By default, population size 10,000 is used."<<endl;
	cout<<setw(20)<<"-mm INPUT"<<"  --  "<<"Multiple merger parameters are defined by a single numerical constant, "<<endl;
	cout<<setw(26)<<" "<<"or a string which speifies the parameter on each branch. "<<endl;
	cout<<setw(26)<<" "<<"The string can be input through command line or a file. "<<endl;
	cout<<setw(26)<<" "<<"By default, Kingman coalescent is used."<<endl;
	cout<<setw(20)<<"-S n1 n2 ..."<<"  --  "<<"Specify the number of samples for each taxon."<<endl;
	cout<<setw(20)<<"-num N"<<"  --  "<<"The number of gene trees will be simulated."<<endl;
	cout<<setw(20)<<"-seed SEED"<<"  --  "<<"User define random SEED"<<endl;
	cout<<setw(20)<<"-mu MU"<<"  --  "<<"User defined constant mutation rate MU. By default mutation rate 0.00005 is used."<<endl;
	cout<<setw(20)<<"-o FILE [option]"<<"  --  "<<"Specify the file name prefix for simulated gene trees. \"GENE_TREE\" by default"<<endl;
	//cout<<"     By default, gene tree branch lengths are in coalescent unit "<<endl;
	cout<<setw(20)<<"-sim_mut_unit"<<"  --  "<<"Convert the simulated gene tree branch lengths to mutation unit."<<endl;
	cout<<setw(20)<<"-sim_num_gener"<<"  --  "<<"Convert the simulated gene tree branch lengths to number of generations."<<endl;
	cout<<setw(20)<<"-sim_num_mut"<<"  --  "<<"Simulate numbers of mutations on each branch of simulated gene trees."<<endl;
	cout<<setw(20)<<"-sim_Si_num"<<"  --  "<<"Generate the file out table, which includes the number of segregating"<<endl;
	cout<<setw(26)<<" "<<"sites and the total branch length of the gene tree in coalescent unit."<<endl;
	cout<<setw(20)<<"-f"<<"  --  "<<"Generate a topology frequency table of a set of input trees or simulated gene trees."<<endl;
	cout<<setw(26)<<" "<<"Frequency table is saved in file freq out by default."<<endl;
	cout<<setw(20)<<"-fF FILE"<<"  --  "<<"The topology frequency table will be saved in the FILE."<<endl;
	cout<<setw(20)<<"-gt FILE"<<"  --  "<<"Specify the FILE of trees to analyse tree topology frequencies."<<endl;
	cout<<setw(20)<<"-mono"<<"  --  "<<"Generate a frequency table of monophyletic, paraphyletic and polyphyletic trees. "<<endl;
	cout<<setw(20)<<"-plot/-dot [option]"<<"  --  "<<"Use LaTEX(-plot) or Dot (-dot) to draw the input (defined by -spcu) network(tree)."<<endl;
	cout<<setw(20)<<"      -branch"<<"  --  "<<"Branch lengths will be labelled in the figure."<<endl;
	cout<<setw(20)<<"-plotF/-dotF FILE"<<"  --  "<<"Generated figure will be saved in FILE."<<endl;			
	cout<<endl;	
}

/*! \brief hybrid-Lambda help file*/
void print_help(){
	cout<<endl;
	cout<<endl;
	cout<<"*****************************************************************"<<endl;
	cout<<"*                      hybrid-Lambda beta 0.4                   *"<<endl;
	cout<<"*                         Author: Joe ZHU                       *"<<endl;
	cout<<"*****************************************************************"<<endl;
	cout<<endl<<endl;
	print_option();
	print_example();
    exit (EXIT_SUCCESS);
}



