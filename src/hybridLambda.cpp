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

// parameters
#include<hybridLambda.hpp>
void HybridLambda::init(){
	this->seed            = (unsigned)(time(0));
	this->simulation_bool = false;
	this->help            = false;
	this->freq_bool       = false;
	this->print_tree      = false;
	this->plot_bool       = false;
	this->log_bool        = false;
	//this->log_NAME        = "LOG";
	this->seg_bool        = false;
	this->read_GENE_trees = false;
	this->read_mt_trees   = false;
	//this->tmrca_NAME      = "tmrcaFILE";
	//this->bl_NAME         = "blFILE";
    this->firstcoal_bool  = false;
    //this->firstcoal_NAME  = "firstcoalFILE" ; 
    this->fst_bool        = false;
}

//HybridLambda::param(){
	//this->init();
//};

HybridLambda::HybridLambda(int argc, char *argv[]){
	this->init();
	int argc_i=1;
	while (argc_i < argc){	
		std::string argv_i(argv[argc_i]);
		if (argv_i=="-h" || argv_i=="-help"){
			help=true;
			print_help();
		}
		
		if (argv_i=="-sp_coal_unit" || argv_i=="-sp_num_gener" || argv_i=="-spcu" || argv_i=="-spng"){
			simulation_bool=true;
		}
		
		if (argv_i=="-seed"){
			read_input_to_param<size_t>(argv[argc_i+1],seed);
			argc_i++;
		}
			
		if (argv_i=="-gt"){
			read_GENE_trees=true;
			gt_file_name=argv[argc_i+1];
			argc_i++;
		}

		if (argv_i=="-mt"){/*! read number of mutations site and simulate segregating sites*/
			read_mt_trees=true;
			mt_file_name=argv[argc_i+1];
			argc_i++;

		}

		if (argv_i=="-mm"){
			mm_bool=true;
		}
			
		if (argv_i=="-pop"){
			pop_bool=true;			
		}
		
		if (argv_i=="-tmrca"){
			tmrca_bool=true;
			//argc_i++;
			//if (argc_i < argc){
				//if (argv[argc_i][0]!='-'){
					//tmrca_NAME=argv[argc_i];
					////argc_i++;
				//}
				//else{argc_i--;}
			//}
		}

		if (argv_i=="-firstcoal"){
			firstcoal_bool=true;
			//argc_i++;
			//if (argc_i < argc){
				//if (argv[argc_i][0]!='-'){
					//firstcoal_NAME=argv[argc_i];
					////argc_i++;
				//}
				//else{argc_i--;}
			//}
		}

		if (argv_i=="-bl"){
			bl_bool=true;
			//argc_i++;
			//if (argc_i < argc){
				//if (argv[argc_i][0]!='-'){
					//bl_NAME=argv[argc_i];
					////argc_i++;
				//}
				//else{argc_i--;}
			//}
		}

		if ( argv_i == "-freq"|| argv_i == "-f" || argv_i == "-freq_file"|| argv_i == "-fF" ){
			freq_bool=true;
		}
		
		if ( argv_i == "-plot" || argv_i == "-plot_file" || argv_i == "-plotF" || argv_i == "-dot" || argv_i == "-dot_file" || argv_i == "-dotF"){
			plot_bool=true;
		}
        
        if ( argv_i == "-label" || argv_i == "-branch" ){
            continue;
            }

		if (argv_i=="-seg" || argv_i=="-segD"){
			seg_bool=true;
			//sim_num_mut_bool=true;
		}
				
		if (argv_i=="-print"){
			print_tree=true;
		}

		if (argv_i=="-fst"){
			fst_bool = true;
		}		

		//if (argv_i=="-log"){
			//log_bool=true;
			//argc_i++;
			//if (argc_i < argc){
				//if (argv[argc_i][0]!='-'){
					//log_NAME=argv[argc_i];
					////argc_i++;
				//}
				//else{argc_i--;}
			//}
		//}
        argc_i++;
	
	}
		
	//srand(seed);	// initialize gnu seed
	MTRand_closed mt;
	mt.seed(seed);		// initialize mt seed
}


void HybridLambda::extract_tmrca(){
    if ( !this->tmrca_bool ){ return; }
    
    this->extract_file_name = this->prefix + "tmrca";
    remove ( this->extract_file_name.c_str() );
    this->extract_file.open ( this->extract_file_name.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
    for ( size_t i = 0; i < gt_tree_str_s.size(); i++ ){
        Net gt(gt_tree_str_s[i]);
        this->extract_file << gt.Net_nodes.back().absolute_time << endl;
    }
    this->extract_file.close();
    std::clog << "TMRCA file is saved at: "<< extract_file_name << "\n";
}


void HybridLambda::extract_bl(){
    if ( !this->bl_bool ){ return; }

    this->extract_file_name = this->prefix + "BL";
    remove ( this->extract_file_name.c_str() );
    this->extract_file.open ( this->extract_file_name.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
    for ( size_t i=0; i < gt_tree_str_s.size(); i++ ){
        Net gt(gt_tree_str_s[i]);
        double totalbl=0;
        for (size_t node_i = 0 ; node_i < gt.Net_nodes.size(); node_i++){
            totalbl = totalbl + gt.Net_nodes[node_i].brchlen1;
        }
        this->extract_file << totalbl << endl;
    }
    this->extract_file.close();
    std::clog << "Total branch length file is saved at: "<< extract_file_name << "\n";
}


void HybridLambda::extract_firstcoal(){
    if ( !this->firstcoal_bool ){ return; }
    
    this->extract_file_name = this->prefix+"fistcoal";
    remove ( this->extract_file_name.c_str() );
    this->extract_file.open ( this->extract_file_name.c_str(), std::ios::out | std::ios::app | std::ios::binary); 
    for ( size_t i = 0; i < gt_tree_str_s.size(); i++ ){
        Net gt(gt_tree_str_s[i]);
        size_t first_coal_index_dummy = gt.first_coal_index();
        extract_file << gt.Net_nodes[first_coal_index_dummy].absolute_time << "\t" << gt.Net_nodes[first_coal_index_dummy].clade << endl;
    }
    this->extract_file.close();
    std::clog << "First Coalescent event is saved at: "<< extract_file_name << "\n";

}


/*! \brief  sim_n_gt constructor */
//sim_n_gt::sim_n_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string, vector < int > sample_size,double mutation_rate,int num_sim_gt,bool sim_mut_unit_bool, bool sim_num_gener_bool,bool sim_num_mut_bool,bool mono_bool){
//sim_n_gt::sim_n_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string, vector < int > sample_size,double mutation_rate,int num_sim_gt,action_board my_action){
void HybridLambda::HybridLambda_core(sim::param sim_param,action_board my_action){
	dout<<" Start simulating "<<sim_param.num_sim_gt <<" gene trees -- begin of sim_n_gt::sim_n_gt(sim::param sim_param,action_board my_action)"<<endl;
	
    string para_string=sim_param.para_string;
    
	//string gene_tree_file=my_action.gene_tree_file;
	string sp_string_coal_unit=sim_param.sp_string_coal_unit;
	string sp_string_pop_size=sim_param.sp_string_pop_size;
	
	vector < int > sample_size=sim_param.sample_size;
	double mutation_rate=sim_param.mutation_rate;
	int num_sim_gt=sim_param.num_sim_gt;
	
	string gene_tree_file_coal_unit = this->prefix + "_coal_unit";
	string gene_tree_file_mut_unit  = this->prefix + "_mut_unit";
	string gene_tree_file_num_gener = this->prefix + "_num_gener";
	string gene_tree_file_num_mut   = this->prefix + "_num_mut";
    
	remove(gene_tree_file_coal_unit.c_str());
	remove(gene_tree_file_mut_unit.c_str());
	remove(gene_tree_file_num_gener.c_str());
	remove(gene_tree_file_num_mut.c_str());
	
   	sim_gt_file_coal_unit.open ( gene_tree_file_coal_unit.c_str(), ios::out | ios::app | ios::binary); 
    sim_gt_file_mut_unit.open  ( gene_tree_file_mut_unit.c_str(),  ios::out | ios::app | ios::binary); 
    sim_gt_file_num_gener.open ( gene_tree_file_num_gener.c_str(), ios::out | ios::app | ios::binary); 
    sim_gt_file_num_mut.open   ( gene_tree_file_num_mut.c_str(),   ios::out | ios::app | ios::binary); 

	//if (my_action.sim_mut_unit_bool){
		//sim_gt_file_mut_unit.open (gene_tree_file_mut_unit.c_str(), ios::out | ios::app | ios::binary); 
	//}
	//if (my_action.sim_num_gener_bool){
		//sim_gt_file_num_gener.open (gene_tree_file_num_gener.c_str(), ios::out | ios::app | ios::binary); 
	//}
	//if (my_action.sim_num_mut_bool){
		//sim_gt_file_num_mut.open (gene_tree_file_num_mut.c_str(), ios::out | ios::app | ios::binary); 
	//}

	if (my_action.Si_num_bool){
		int total_lineage=0;
		for (size_t i=0; i<sim_param.sample_size.size();i++){
			total_lineage=total_lineage+sim_param.sample_size[i];
		}
		outtable_header(total_lineage);
	}
	for (int i=0;i<num_sim_gt;i++){
		//sim_one_gt sim_gt_string(sp_string_coal_unit, sp_string_pop_size, para_string, sample_size,  mutation_rate,my_action);
		sim_one_gt sim_gt_string(sim_param,my_action);
		gt_tree_str_s.push_back(sim_gt_string.gt_string_coal_unit);
		if (my_action.sim_num_mut_bool){
			mt_tree_str_s.push_back(sim_gt_string.gt_string_mut_num);
		}
		
		if (my_action.mono_bool){
			if (i==0){
				tax_name=sim_gt_string.tax_name;
				monophyly=sim_gt_string.monophyly;	
			}
			else{
				for (unsigned int mono_i=0;mono_i<monophyly.size();mono_i++){
					monophyly[mono_i]=monophyly[mono_i]+sim_gt_string.monophyly[mono_i];
				}		
			}
		}
		
		dout << sim_gt_string.gt_string_coal_unit << endl;
		sim_gt_file_coal_unit<<sim_gt_string.gt_string_coal_unit <<"\n";
		        
		if ( my_action.sim_mut_unit_bool ){ sim_gt_file_mut_unit  << sim_gt_string.gt_string_mut_unit  << "\n"; }
		if ( my_action.sim_num_gener_bool){ sim_gt_file_num_gener << sim_gt_string.gt_string_gener_num << "\n"; }
		if ( my_action.sim_num_mut_bool  ){ sim_gt_file_num_mut   << sim_gt_string.gt_string_mut_num   << "\n"; }
	}
	
	if (my_action.mono_bool){
		for (unsigned int mono_i=0;mono_i<monophyly.size();mono_i++){
			monophyly[mono_i]=monophyly[mono_i]/num_sim_gt;
		}
	}
    
    sim_gt_file_coal_unit.close();
    sim_gt_file_mut_unit.close();
    sim_gt_file_num_gener.close();
    sim_gt_file_num_mut.close();
	//if (my_action.sim_mut_unit_bool){
		//sim_gt_file_mut_unit.close();
	//}
	//if (my_action.sim_num_gener_bool){
		//sim_gt_file_num_gener.close();
	//}
	//if (my_action.sim_num_mut_bool){
		//sim_gt_file_num_mut.close();
	//}
	dout<<"end of sim_n_gt::sim_n_gt(sim::param sim_param,action_board my_action)"<<endl;

}


/*! \brief hybrid-Lambda help file*/
void HybridLambda::print_help(){
	cout<<endl;
	cout<<endl;
	cout<<"*****************************************************************"<<endl;
	cout<<"*                      hybrid-Lambda beta 0.4                   *"<<endl;
	cout<<"*                         Author: Joe ZHU                       *"<<endl;
	cout<<"*****************************************************************"<<endl;
	cout<<endl<<endl;
	HybridLambda::print_option();
	HybridLambda::print_example();
	exit(1);
}

void HybridLambda::print_example(){
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

void HybridLambda::print_option(){
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



