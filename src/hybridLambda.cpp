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

#include "figure.hpp"

#include<hybridLambda.hpp>
#include <stdlib.h>     /* strtod */

void HybridLambda::init(){
	this->seed            = (unsigned)(time(0));
	this->simulation_bool = false;
	this->prefix = "OUT";
	this->freq_bool       = false;
	this->print_tree_bool = false;
	this->plot_bool       = false;
	this->seg_bool        = false;
	this->read_GENE_trees = false;
	this->read_mt_trees   = false;
    this->firstcoal_bool  = false;
    this->fst_bool        = false;
    
    this->tmrca_bool = false;
    this->bl_bool = false;
    this->firstcoal_bool = false;
    this->argc_i = 1;
    
    this->gt_file_name = "";
    this->mt_file_name = "";
    this->num_sim_gt = 1;
    this->simulation_jobs_ = new action_board();
    this->parameters_ = new SimulationParameters();
}

HybridLambda::~HybridLambda(){ 
    delete simulation_jobs_; 
    delete parameters_; 
}

string HybridLambda::read_input_para(const char *inchar, string in_str){
    if ( is_num(inchar) ) return write_para_into_tree( in_str, strtod(inchar, NULL) );
    else return read_input_line( inchar );
}

void HybridLambda::read_sp_str( string & argv_i ){
    this->simulation_bool=true; 
    if (argv_i == "-sp_num_gener" || argv_i == "-spng"){ this->parameters_->num_gener_bool=true; }
    if (argv_i == "-sp_coal_unit" || argv_i == "-spcu"){ this->parameters_->sp_coal_unit_bool=true;	}
    if ( this->parameters_->sp_coal_unit_bool && this->parameters_->num_gener_bool){
        throw std::invalid_argument("Species tree branch length should only be in Coalescent unit or number of generations. Choose either -sp_num_gener or -sp_coal_unit");
    }
    readNextStringto( this->tmp_input_str , this->argc_i, this->argc_,  this->argv_ );
    this->parameters_->net_str = read_input_line(tmp_input_str.c_str());
}

void HybridLambda::extract_mm_or_pop_param( string & mm_pop_string ){
    readNextStringto( this->tmp_input_str , this->argc_i, this->argc_,  this->argv_ );
    mm_pop_string = this->read_input_para(this->tmp_input_str.c_str(), this->parameters_->net_str);
}


void HybridLambda::parse(){
	while (argc_i < argc_){	
		std::string argv_i(argv_[argc_i]);
		if ( argv_i == "-h" || argv_i == "-help" ){ print_help(); }
		else if ( argv_i =="-gt" ){ readNextStringto( this->gt_file_name , this->argc_i, this->argc_,  this->argv_ ); }
        /*! read number of mutations site and simulate segregating sites*/
		else if ( argv_i =="-mt" ){ readNextStringto( this->mt_file_name , this->argc_i, this->argc_,  this->argv_ ); }
		else if ( argv_i =="-seed" ){ this->seed = this->readNextInput<size_t>(); }
        else if ( argv_i =="-num"){ this->num_sim_gt = this->readNextInput<int>(); }
        else if ( argv_i =="-sp_coal_unit" || argv_i=="-sp_num_gener" || argv_i == "-spcu" || argv_i=="-spng"){ this->read_sp_str(argv_i); }
		else if ( argv_i =="-mm" ){ this->parameters_->mm_bool = true; this->extract_mm_or_pop_param( this->parameters_->para_string ) ; }
        else if ( argv_i =="-pop" ){ this->parameters_->pop_bool = true; this->extract_mm_or_pop_param( this->parameters_->sp_string_pop_size ) ; }        
        else if ( argv_i =="-S" ){ this->read_sample_sizes();	}
        else if ( argv_i =="-mu"){ this->parameters_->mutation_rate = this->readNextInput<double>(); }	
        else if ( argv_i =="-o" ) readNextStringto( this->prefix , this->argc_i, this->argc_,  this->argv_ );
        else if ( argv_i == "-sim_mut_unit"  ){ this->simulation_jobs_->set_sim_mut_unit(); }
        else if ( argv_i == "-sim_num_gener" ){ this->simulation_jobs_->set_sim_num_gener();}
		else if ( argv_i == "-sim_num_mut" ){ this->simulation_jobs_->set_sim_num_mut(); } 
		else if ( argv_i == "-sim_Si_num"    ){ this->simulation_jobs_->set_Si_num(); check_and_remove("out_table");} // work on code for removing out_table
		else if ( argv_i == "-mono"){ this->simulation_jobs_->set_mono() ; }
		else if ( argv_i == "-freq" || argv_i == "-f" ){ this->freq_bool=true; }
        //else if ( argv_i == "-freq_file"|| argv_i == "-fF" ){ this->freq_bool=true;  this->argc_i++;}
		else if ( argv_i == "-plot" || argv_i == "-dot" ){ this->plot_bool = true; }
        //else if ( argv_i == "-plot_file" || argv_i == "-plotF" || argv_i == "-dot_file" || argv_i == "-dotF" ){ this->plot_bool=true; this->argc_i++; }
        else if ( argv_i == "-label" || argv_i == "-branch" ){ argc_i++; continue; }        
		else if ( argv_i == "-tmrca"     ){ this->tmrca_bool = true;     }
		else if ( argv_i == "-firstcoal" ){ this->firstcoal_bool = true; }
		else if ( argv_i == "-bl"   ){ this->bl_bool = true; }       
		else if ( argv_i == "-seg"  ){ this->seg_bool = true;	}
        //else if ( argv_i == "-segD" ){ this->seg_bool=true; readNextStringto( this->seg_dir_name , this->argc_i, this->argc_,  this->argv_ ); }        
		else if ( argv_i == "-print" ){ this->print_tree_bool = true; }
		else if ( argv_i == "-fst"   ){ this->fst_bool = true; }
		else { throw std::invalid_argument ( "Unknown flag:" + argv_i); }
        //else { cout <<"  need to change this !!!" << argv_i<<endl; argc_i++;continue; } // need to change this !!!
        argc_i++;	
	}
		
    this->finalize();
}

void HybridLambda::read_sample_sizes(){
    this->parameters_->samples_bool = true;
    while ( argc_i < argc_ ){
        try {
            int pop_num = this->readNextInput<int>();
            this->parameters_->sample_size.push_back(pop_num);
        } catch (std::invalid_argument e) {
            --argc_i;
            break;
        }
    }
}


void HybridLambda::finalize(){
    if ( this->gt_file_name.size() > 0 ) this->read_input_lines( this->gt_file_name.c_str(), this->gt_tree_str_s);
    if ( this->mt_file_name.size() > 0 ) this->read_input_lines( this->mt_file_name.c_str(), this->mt_tree_str_s);
    if ( this->simulation_bool ) this->parameters()->finalize( );    
    if ( this->seg_bool ) {
        this->simulation_jobs_->set_sim_num_mut();
        this->seg_dir_name = this->prefix + "seg_sites" ;
        }
    if ( this->print_tree_bool ) this->print();
    if ( this->plot_bool ){
        Figure figure_para ( this->argc_, this->argv_ );
        figure_para.figure_file_prefix = this->prefix;
        figure_para.finalize();
        figure_para.plot( this->parameters()->net_str );
        exit(EXIT_SUCCESS);
        }
}

void HybridLambda::print(){
    Net net( this->parameters_->net_str );
    net.print_all_node();
    exit(EXIT_SUCCESS);
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
void HybridLambda::HybridLambda_core( ){
    if ( !simulation_bool ){ return; }
   	//srand(seed);	// initialize gnu seed
    MTRand_closed mt;
	mt.seed( this->seed );		// initialize mt seed
    clog << "Random seed: " << this->seed <<endl;
	dout<<" Start simulating "<< this->num_sim_gt <<" gene trees -- begin of sim_n_gt::sim_n_gt(sim::param sim_param,action_board my_action)"<<endl;
		
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

	if ( this->simulation_jobs_->Si_num_bool ){
		int total_lineage=0;
		for ( size_t i = 0; i < this->parameters_->sample_size.size();i++){
			total_lineage=total_lineage+ this->parameters_->sample_size[i];
		}
		outtable_header(total_lineage);
	}
	for ( int i=0; i < this->num_sim_gt; i++ ){
		sim_one_gt sim_gt_string( this->parameters_, this->simulation_jobs_ );
		gt_tree_str_s.push_back(sim_gt_string.gt_string_coal_unit);
		if ( this->simulation_jobs_->sim_num_mut_bool){
			mt_tree_str_s.push_back(sim_gt_string.gt_string_mut_num);
		}
		
		if ( this->simulation_jobs_->mono_bool){
			if ( i == 0 ){
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
		        
		if (  this->simulation_jobs_->sim_mut_unit() ){ sim_gt_file_mut_unit  << sim_gt_string.gt_string_mut_unit  << "\n"; }
		if (  this->simulation_jobs_->sim_num_gener()){ sim_gt_file_num_gener << sim_gt_string.gt_string_gener_num << "\n"; }
		if (  this->simulation_jobs_->sim_num_mut()  ){ sim_gt_file_num_mut   << sim_gt_string.gt_string_mut_num   << "\n"; }
	}
	
	if ( this->simulation_jobs_->mono()){
		for ( size_t mono_i = 0; mono_i < monophyly.size(); mono_i++ ){
			monophyly[mono_i] = monophyly[mono_i] / num_sim_gt;
		}
	}
    
    sim_gt_file_coal_unit.close();
    sim_gt_file_mut_unit.close();
    sim_gt_file_num_gener.close();
    sim_gt_file_num_mut.close();
    
    std::clog << "Produced gene tree files: \n";	
	std::clog << gene_tree_file_coal_unit << "\n";				
    if (  this->simulation_jobs_->sim_mut_unit_bool  ) std::clog << gene_tree_file_mut_unit  << "\n"; 
    else remove (gene_tree_file_mut_unit.c_str()) ;
    if (  this->simulation_jobs_->sim_num_gener_bool ) std::clog << gene_tree_file_num_gener << "\n"; 
    else remove (gene_tree_file_num_gener.c_str()) ;
    if (  this->simulation_jobs_->sim_num_mut_bool   ) std::clog << gene_tree_file_num_mut   << "\n"; 
    else remove (gene_tree_file_num_mut.c_str()) ;
	dout<<"end of sim_n_gt::sim_n_gt(sim::param sim_param,action_board my_action)"<<endl;
}


string HybridLambda::read_input_line(const char *inchar){
	string out_str;
    ifstream in_file( inchar );
	if (in_file.good())	getline ( in_file, out_str); 
	else{
		string dummy_str(inchar);
		if (dummy_str.find('(')!=string::npos && dummy_str.find(')')!=string::npos) out_str=dummy_str;
		else  throw std::invalid_argument("Invalid input file. " + string (inchar) );
	}
	in_file.close();			
    return 	out_str;
}

void HybridLambda::read_input_lines(const char inchar[], vector <string> & out_vec){
	ifstream in_file( inchar );
	string out_str;
	if ( in_file.good() ){
		getline ( in_file, out_str );
		while ( out_str.size() > 0 ){   
			out_vec.push_back( out_str );
			getline ( in_file, out_str );
		}
	}	
	else{
		string dummy_str(inchar);
		if (dummy_str.find('(')!=string::npos && dummy_str.find(')')!=string::npos){
			out_str=dummy_str;
			out_vec.push_back(out_str);
		}else{
			throw std::invalid_argument("Invalid input file. " + string (inchar) );
		}
	}
	in_file.close();
}


bool HybridLambda::is_num(const char *inchar){
	bool is_num_return=true;
	string in_str(inchar);
	for (size_t i=0;i<in_str.size();i++){
		if (isalpha(in_str[i]) && in_str[i]!='e'){
			is_num_return=false;
			break;
		}
	}
	return is_num_return;
}

/*! \brief remove old segregating sites data, and generate new ones */
void HybridLambda::create_site_data_dir(){
    if ( !this->seg_bool ){ return; }
	string rm_commond="rm -rf " + this->seg_dir_name;
	int sys = system( rm_commond.c_str() );
	string mkdir_commond="mkdir "+ this->seg_dir_name;
	sys=system( mkdir_commond.c_str() );
	for ( size_t i = 0; i < mt_tree_str_s.size(); i++ ){
		create_new_site_data( mt_tree_str_s[i], i+1 );
	}
    clog << "Segregating site data saved at: "<< this->seg_dir_name <<"\n"; 
}

/*! \brief Generate segrateing site data */
void HybridLambda::create_new_site_data( string &gt_string_mut_num, int site_i ){
	Net mt_tree( gt_string_mut_num );
    ostringstream site_i_str;
	site_i_str << site_i;
	string sitefile_name = seg_dir_name + "/site" + site_i_str.str();
	extract_file.open (sitefile_name.c_str()); 
	
	int total_mut = 0;
	for ( size_t node_i = 0; node_i < mt_tree.Net_nodes.size(); node_i++ ){
		total_mut = total_mut + mt_tree.Net_nodes[node_i].brchlen1;
	}
	for ( size_t tip_i = 0; tip_i < mt_tree.tip_name.size(); tip_i++ ){
		extract_file << mt_tree.tip_name[tip_i] << " ";
		for ( size_t node_i = 0; node_i < mt_tree.Net_nodes.size(); node_i++ ){
			if ( mt_tree.Net_nodes[node_i].brchlen1 > 0 ){
				for ( int num_repeat = 0; num_repeat < mt_tree.Net_nodes[node_i].brchlen1; num_repeat++ ){				
					extract_file << mt_tree.descndnt2[node_i][tip_i] ;
				}
			}
		}
		extract_file<<"\n";
	}
	extract_file.close();	
}
