/*
 * hybrid-Lambda is used to simulate gene trees given species network under
 * coalescent process.
 *
 * Copyright (C) 2010 -- 2015 Sha (Joe) Zhu
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

#include "hybridLambda.hpp"
#include "plot/figure.hpp"
#include "freq/freq.hpp"
#include <stdlib.h>     /* strtod */
#include "mersenne_twister.hpp"


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
    this->fst_bool        = false;
    this->tmrca_bool      = false;
    this->bl_bool         = false;
    this->firstcoal_bool  = false;
    this->print_help_bool = false;
    this->argc_i = 1;

    this->gt_file_name = "";
    this->mt_file_name = "";
    this->num_sim_gt = 1;
    this->simulation_jobs_ = new action_board();
    this->parameters_      = new SimulationParameters();
}

HybridLambda::~HybridLambda(){
    delete simulation_jobs_;
    delete parameters_;
}

string HybridLambda::read_input_para( const char *inchar, string in_str ){
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
    if ( argc_ == 1 ){ this->print_help_bool = true; }

    while (argc_i < argc_){
        std::string argv_i(argv_[argc_i]);
        if ( argv_i == "-h" || argv_i == "-help" ){ this->print_help_bool = true; }
        else if ( argv_i =="-gt" ){ readNextStringto( this->gt_file_name , this->argc_i, this->argc_,  this->argv_ ); }
        /*! read number of mutations site and simulate segregating sites*/
        else if ( argv_i =="-mt" ){ readNextStringto( this->mt_file_name , this->argc_i, this->argc_,  this->argv_ ); }
        else if ( argv_i =="-seed" ){ this->seed = this->readNextInput<size_t>(); }
        else if ( argv_i =="-num" ){ this->num_sim_gt = this->readNextInput<int>(); }
        else if ( argv_i =="-sp_coal_unit" || argv_i=="-sp_num_gener" || argv_i == "-spcu" || argv_i=="-spng"){ this->read_sp_str(argv_i); }
        else if ( argv_i =="-mm" ){ this->parameters_->mm_bool = true; this->extract_mm_or_pop_param( this->parameters_->para_string ) ; }
        else if ( argv_i =="-pop" ){ this->parameters_->pop_bool = true; this->extract_mm_or_pop_param( this->parameters_->sp_string_pop_size ) ; }
        else if ( argv_i =="-S" ){ this->read_sample_sizes();	}
        //else if ( argv_i =="-mu"){ this->parameters_->mutation_rate = this->readNextInput<double>(); }
        else if ( argv_i =="-mu"){ ++argc_i; this->parameters_->mutation_rate = strtod(argv_[argc_i], NULL); }
        else if ( argv_i =="-o" ) readNextStringto( this->prefix , this->argc_i, this->argc_,  this->argv_ );
        else if ( argv_i == "-sim_mut_unit"  ){ this->simulation_jobs_->set_sim_mut_unit(); }
        else if ( argv_i == "-sim_num_gener" ){ this->simulation_jobs_->set_sim_num_gener();}
        else if ( argv_i == "-sim_num_mut" ){ this->simulation_jobs_->set_sim_num_mut(); }
        else if ( argv_i == "-sim_Si_num"    ){ this->simulation_jobs_->set_Si_num(); } // work on code for removing out_table
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
    if ( this->print_help_bool ){
        this->print_help();
        return;
    }

    if ( this->simulation_bool ) {
        if ( this->fst_bool ){
            this->fst_file_name = this->prefix + "_fst";
            remove ( fst_file_name.c_str() );
        }
        if ( this->gt_file_name.size() > 0 ) clog << "WARNING: \"-gt\" option is ignored" << endl;
        if ( this->mt_file_name.size() > 0 ) clog << "WARNING: \"-mt\" option is ignored" << endl;

        this->parameters()->finalize( );

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
            //exit(EXIT_SUCCESS);
            return;
        }
    } else {
        if ( this->gt_file_name.size() > 0 ) {
            this->read_input_lines( this->gt_file_name.c_str(), this->gt_tree_str_s);
        } else if ( this->mt_file_name.size() > 0 ) {
            this->read_input_lines( this->mt_file_name.c_str(), this->mt_tree_str_s);
            this->seg_dir_name = this->prefix + "seg_sites" ;  // Initialize segregating site data directory
        } else {
            delete simulation_jobs_;
            delete parameters_;
            throw std::invalid_argument( "No input was provided!" );
        }
    }
}

void HybridLambda::print(){
    Net net( this->parameters_->net_str );
    net.print_all_node();
    return;
    //exit(EXIT_SUCCESS);
}

void HybridLambda::extract_tmrca(){
    if ( !this->tmrca_bool ){ return; }

    this->extract_file_name = this->prefix + "tmrca";
    remove ( this->extract_file_name.c_str() );
    this->extract_file.open ( this->extract_file_name.c_str(), std::ios::out | std::ios::app | std::ios::binary);
    for ( size_t i = 0; i < gt_tree_str_s.size(); i++ ){
        Tree gt(gt_tree_str_s[i]);
        this->extract_file << gt.NodeContainer.back().height() << endl;
    }
    this->extract_file.close();
    std::clog << "TMRCA file is saved at: "<< extract_file_name << "\n";
}


void HybridLambda::extract_bl(){
    if ( !this->bl_bool ){ return; }

    this->extract_file_name = this->prefix + "bl";
    remove ( this->extract_file_name.c_str() );
    this->extract_file.open ( this->extract_file_name.c_str(), std::ios::out | std::ios::app | std::ios::binary);
    for ( size_t i=0; i < gt_tree_str_s.size(); i++ ){
        Tree gt(gt_tree_str_s[i]);
        double totalbl = 0;
        vector <double> bl_of_x_descndnt(gt.tax_name.size() - 1, 0.0);
        for (size_t node_i = 0 ; node_i < gt.NodeContainer.size(); node_i++){
            totalbl += gt.NodeContainer[node_i].brchlen1();
            bl_of_x_descndnt[gt.NodeContainer[node_i].num_descndnt - 1] += gt.NodeContainer[node_i].brchlen1();
        }

        this->extract_file << totalbl << "\t";
        for ( auto const &value: bl_of_x_descndnt){
            this->extract_file << value << "\t";
        }

        this->extract_file << endl;
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
        Tree gt(gt_tree_str_s[i]);
        size_t first_coal_index_dummy = gt.first_coal_index();
        extract_file << gt.NodeContainer[first_coal_index_dummy].height() << "\t" << gt.NodeContainer[first_coal_index_dummy].clade << endl;
    }
    this->extract_file.close();
    std::clog << "First Coalescent event is saved at: "<< extract_file_name << "\n";

}

void HybridLambda::extract_frequency(){
    if ( !this->freq_bool ) return;
    Frequency freq_para;
    freq_para.freq_out_filename = this->prefix + "_frequencies";

    ifstream tmp_file( freq_para.freq_out_filename.c_str() );
	if ( tmp_file.good() ) 	{  remove(freq_para.freq_out_filename.c_str()); 	}

    freq_para.compute_gt_frequencies( this->gt_tree_str_s );
}


bool HybridLambda::mono_fst_not_feasiable( string flag ){ // flag is either -mono or -fst
    if ( this->parameters()->sample_size.size() != 2 || this->parameters()->is_Net ) {
        std::clog << "ERROR: \"" + flag + "\" flag only applies to species tree of two population." << endl;
        return true;
    }
    else return false;
}

void HybridLambda::extract_mono() {
    if ( !this->simulation_jobs()->mono() ) { return ;}
    if ( this->mono_fst_not_feasiable ("-mono") ) { return; }

    cout << "   A mono     B mono Recip mono     A para     B para  Polyphyly" << endl;
    for ( size_t mono_i = 0; mono_i < this->monophyly.size(); mono_i++){
        cout << setw(9) << this->monophyly[mono_i] / this->num_sim_gt << "  ";
    }
    cout << endl;
}

/*! \brief  sim_n_gt constructor */
void HybridLambda::HybridLambda_core( ){

    if ( !simulation_bool ){ return; }
    MersenneTwister rg(this->seed);
    clog << "Random seed: " << this->seed <<endl;
    dout <<" Start simulating "<< this->num_sim_gt <<" gene trees -- begin of sim_n_gt::sim_n_gt(sim::param sim_param,action_board my_action)"<<endl;

    string gene_tree_file_coal_unit = this->prefix + "_coal_unit";
    string gene_tree_file_mut_unit  = this->prefix + "_mut_unit";
    string gene_tree_file_num_gener = this->prefix + "_num_gener";
    string gene_tree_file_num_mut   = this->prefix + "_num_mut";
    string gene_tree_file_Si_table  = this->prefix + "_Si_table";

    remove(gene_tree_file_coal_unit.c_str());
    remove(gene_tree_file_mut_unit.c_str());
    remove(gene_tree_file_num_gener.c_str());
    remove(gene_tree_file_num_mut.c_str());
    remove(gene_tree_file_Si_table.c_str());

    sim_gt_file_coal_unit.open ( gene_tree_file_coal_unit.c_str(), ios::out | ios::app | ios::binary);
    sim_gt_file_mut_unit.open  ( gene_tree_file_mut_unit.c_str(),  ios::out | ios::app | ios::binary);
    sim_gt_file_num_gener.open ( gene_tree_file_num_gener.c_str(), ios::out | ios::app | ios::binary);
    sim_gt_file_num_mut.open   ( gene_tree_file_num_mut.c_str(),   ios::out | ios::app | ios::binary);
    extract_file.open          ( gene_tree_file_Si_table.c_str(),  ios::out | ios::app | ios::binary);

    this->outtable_header(extract_file);

    for ( int i=0; i < this->num_sim_gt; i++ ){
        //this->parameters_->my_Net->print_all_node();
        simTree sim_gt_string( this->parameters_, this->simulation_jobs_ , this->extract_file, &rg);
        gt_tree_str_s.push_back(sim_gt_string.gt_string_coal_unit);
        if ( this->simulation_jobs_->sim_num_mut_bool ) mt_tree_str_s.push_back(sim_gt_string.gt_string_mut_num);

        if ( this->simulation_jobs_->mono_bool){
            if ( i == 0 ){
                monophyly = sim_gt_string.monophyly;
            } else{
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

    sim_gt_file_coal_unit.close();
    sim_gt_file_mut_unit.close();
    sim_gt_file_num_gener.close();
    sim_gt_file_num_mut.close();
    extract_file.close();

    std::clog << "Produced gene tree files: \n";
    std::clog << gene_tree_file_coal_unit << "\n";

    if (  this->simulation_jobs_->sim_mut_unit_bool  ) std::clog << gene_tree_file_mut_unit  << "\n";
    else remove (gene_tree_file_mut_unit.c_str()) ;

    if (  this->simulation_jobs_->sim_num_gener_bool ) std::clog << gene_tree_file_num_gener << "\n";
    else remove (gene_tree_file_num_gener.c_str()) ;

    if (  this->simulation_jobs_->sim_num_mut_bool   ) std::clog << gene_tree_file_num_mut   << "\n";
    else remove (gene_tree_file_num_mut.c_str()) ;

    if (  this->simulation_jobs_->Si_num_bool   ) std::clog << gene_tree_file_Si_table  << "\n";
    else remove (gene_tree_file_Si_table.c_str()) ;

    dout<<"end of sim_n_gt::sim_n_gt(sim::param sim_param,action_board my_action)"<<endl;

    this->extract_mono();
    rg.clearFastFunc(); // Empty memory was allocated for fastfunc
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

/*! \brief Printing out the header for out_table*/
void HybridLambda::outtable_header( std::ofstream &output ){
   	if ( !this->simulation_jobs_->Si_num_bool ){return;}

	output <<"t_total\tt_MRCA\tS_total";
	for ( int sii = 0; sii < this->parameters_->total_num_lineage-1; sii++ ) output<<"\tS_"<<sii+1;
	output << endl;
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
    } else {
        string dummy_str(inchar);
        if (dummy_str.find('(')!=string::npos && dummy_str.find(')')!=string::npos){
            out_str=dummy_str;
            out_vec.push_back(out_str);
        } else{
            throw std::invalid_argument("Invalid input file. " + string (inchar) );
        }
    }
    in_file.close();
}


bool HybridLambda::is_num(const char *inchar){
    bool is_num_return=true;
    string in_str(inchar);
    for (size_t i=0;i<in_str.size();i++){
        if ( isalpha(in_str[i]) && in_str[i]!='e'){
            is_num_return=false;
            break;
        }
    }
    return is_num_return;
}


/*! \brief remove old segregating sites data, and generate new ones */
void HybridLambda::create_site_data_dir(){
    if ( !this->seg_bool ){ return; }
    cout << this->seg_dir_name <<endl;
	string rm_commond="rm -rf " + this->seg_dir_name;
	int sys = system( rm_commond.c_str() );
	string mkdir_commond="mkdir "+ this->seg_dir_name;
	sys=system( mkdir_commond.c_str() );
	for ( size_t i = 0; i < mt_tree_str_s.size(); i++ ){
		create_new_site_data( mt_tree_str_s[i], i+1 );
	}
    clog << "Segregating site data saved at: "<< this->seg_dir_name << "\n";
}

/*! \brief Generate segrateing site data */
void HybridLambda::create_new_site_data( string &gt_string_mut_num, int site_i ){
    Tree mt_tree( gt_string_mut_num );
    string sitefile_name = seg_dir_name + "/site" + std::to_string(static_cast<long long>(site_i));
    extract_file.open (sitefile_name.c_str());

    int total_mut = 0;
    this->haplotypes.clear();
    for ( size_t node_i = 0; node_i < mt_tree.NodeContainer.size(); node_i++ ){
        for ( size_t num_mut = 0 ; num_mut < mt_tree.NodeContainer[node_i].brchlen1(); num_mut++){
            haplotypes.push_back ( mt_tree.samples_below[node_i] );}
        total_mut += mt_tree.NodeContainer[node_i].brchlen1();
    }

    for ( size_t tip_i = 0; tip_i < mt_tree.tip_name.size(); tip_i++ ){
        extract_file << mt_tree.tip_name[tip_i] << "\t";
        for ( size_t i = 0 ; i < haplotypes.size(); i++){
            extract_file << haplotypes[i][tip_i] ;}
        extract_file<<"\n";
    }
    extract_file.close();

    if ( this->fst_bool ) {
        extract_file.open ( fst_file_name.c_str(), std::ios::out | std::ios::app | std::ios::binary );
        extract_file << compute_fst ( this->haplotypes, parameters()->sample_size ) << "\n";
        extract_file.close();
    }
    //cout << "fst: " << compute_fst ( this->haplotypes, parameters()->sample_size ) << endl;
}

void HybridLambda::print_example(){
    cout << "Examples:"
         << endl
         << endl;
    cout << "hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 3 -seed 2 -o example1" << endl;
    cout << "hybrid-Lambda -spcu trees/4_tax_sp_nt1_para -o example2 -num 2 -mu 0.00003 -sim_mut_unit -sim_num_mut" << endl;
    cout << "hybrid-Lambda -spcu '((1:1,2:1):1,3:2);' -num 100 -pop 25000 -sim_num_gener" << endl;
    cout << "hybrid-Lambda -spng '(A:50000,B:50000)r;' -pop '(A:50000,B:50000)r:40000;'" << endl;
    cout << "hybrid-Lambda -spcu '((((A:1.1,B:1.1):2.1,a:2.2):1.1,13D:.2):.3,4:.3);' -S 2 4 3 6 5" << endl;
    cout << "hybrid-Lambda -spcu '(A:1,B:1)r;' -mm '(A:1.9,B:.2)r:2;' -S 3 4" << endl;
    cout << "hybrid-Lambda -spcu trees/7_tax_sp_nt1_para -dot -branch" << endl;
    cout << "hybrid-Lambda -spcu trees/4_tax_sp1.tre -num 1000 -o GENE -f" << endl;
    cout << "hybrid-Lambda -gt GENE_coal_unit -f " << endl;
    cout << "hybrid-Lambda -mt GENE_num_mut -seg " << endl;
    cout << "hybrid-Lambda -spcu '(A:5,B:5)r;' -mono -num 100 -mm .1 -S 4 4" << endl;
    cout << endl;
}

void HybridLambda::print_option(){
    cout << endl
         << "hybrid-Lambda " << VERSION
         << endl
         << endl;
    cout << "Usage:"
         << endl;
    cout << setw(20) << "-h or -help"         << "  --  " << "Help. List the following content." << endl;
    cout << setw(20) << "-spcu STR"           << "  --  " << "Input the species network/tree string through command line or a file." << endl;
    cout << setw(26) << " "                               << "Branch lengths of the INPUT are in coalescent unit." << endl;
    cout << setw(20) << "-spng STR"           << "  --  " << "Input the species network/tree string through command line or a file. " << endl;
    cout << setw(26) << " "                               << "Branch lengths of the INPUT are in number of generation." << endl;
    cout << setw(20) << "-pop STR/FLT"        << "  --  " << "Population sizes are defined by a single numerical constant, " << endl;
    cout << setw(26) << " "                               << "or a string which specifies the population size on each branch. " << endl;
    cout << setw(26) << " "                               << "The string can be input through command line or a file. " << endl;
    cout << setw(26) << " "                               << "By default, population size 10,000 is used." << endl;
    cout << setw(20) << "-mm STR/FLT"         << "  --  " << "Multiple merger parameters are defined by a single numerical constant, " << endl;
    cout << setw(26) << " "                               << "or a string which speifies the parameter on each branch. " << endl;
    cout << setw(26) << " "                               << "The string can be input through command line or a file. " << endl;
    cout << setw(26) << " "                               << "By default, Kingman coalescent is used." << endl;
    cout << setw(20) << "-S INT INT ..."      << "  --  " << "Specify the number of samples for each taxon." << endl;
    cout << setw(20) << "-num INT"            << "  --  " << "The number of gene trees will be simulated." << endl;
    cout << setw(20) << "-seed INT"           << "  --  " << "User define random SEED" << endl;
    cout << setw(20) << "-mu FLT"             << "  --  " << "User defined constant mutation rate per locus. By default mutation rate 0.00005 is used." << endl;
    cout << setw(20) << "-o STR [option]"     << "  --  " << "Specify the file name prefix for simulated gene trees. Prefix is set as \"OUT\" by default." << endl;
    cout << setw(26) << " "                               << "When options are not specified, only output trees with branch lengths are in coalescent unit." << endl;
    cout << setw(20) << "[-sim_mut_unit]"     << "  --  " << "Convert the simulated gene tree branch lengths to mutation unit." << endl;
    cout << setw(20) << "[-sim_num_gener]"    << "  --  " << "Convert the simulated gene tree branch lengths to number of generations." << endl;
    cout << setw(20) << "[-sim_num_mut]"      << "  --  " << "Simulate numbers of mutations on each branch of simulated gene trees." << endl;
    cout << setw(20) << "[-sim_Si_num]"       << "  --  " << "Generate a table, which includes the number of segregating sites" << endl;
    cout << setw(26) << " "                               << "and the total branch length of the gene tree, as well as the TMRCA." << endl;
    cout << setw(20) << "-f"                  << "  --  " << "Generate a topology frequency table of a set of input trees or simulated gene trees." << endl;
    //cout << setw(26) << " " << "Frequency table is saved in file freq out by default." << endl;
    cout << setw(20) << "-gt STR"             << "  --  " << "Specify the FILE NAME of trees to analyse tree topology frequencies." << endl;
    cout << setw(20) << "-seg"                << "  --  " << "Generate segregating site data." << endl;
    cout << setw(20) << "-mt STR"             << "  --  " << "Specify the FILE NAME of trees to generate segregating site data." << endl;
    cout << setw(26) << " "                               << "Tree branch lengths indicate number of mutations on the branch." << endl;
    cout << setw(20) << "-mono"               << "  --  " << "Generate a frequency table of monophyletic, paraphyletic and polyphyletic trees. " << endl;
    cout << setw(20) << "-plot/-dot [option]" << "  --  " << "Use LaTEX(-plot) or Dot (-dot) to draw the input (defined by -spcu) network(tree)." << endl;
    cout << setw(20) << "      -branch"       << "  --  " << "Branch lengths will be labelled in the figure." << endl;
    //cout << setw(20) << "-plotF/-dotF FILE" << "  --  " << "Generated figure will be saved in FILE." << endl;
    cout << endl;
}



/*!
 * Assume two populations A and B have been isolated until time tau in the past as measured from the present.
 * Assume also that the same coalescent process is operating in populations A and B.
 * Let TW denote the time until coalescence for two lines when drawn from the same population,
 * and Tb when drawn from different populations.
 * Let lambdaA denote the coalescence rate for two lines in population A, and
 * lambdaAB for the common ancestral population AB.
 * For the Beta(2 − alpha, alpha)-coalescent, lambdaA = 1, for the point-mass process lambdaA = psi^2. One now obtains
 * ETw exptected value of Tw
 * ETw = (1 - exp(-lambdaA * tau) * lambdaA^{-1} + exp(-lambdaA * tau) * (tau + lambdaAB^{-1})
 *
 */

double lambda( double alpha ){
    return exp(log(binomial_coefficient((double)2,(double)2))+log(Beta(2-alpha,2-2+alpha)) - log(Beta(2.0-alpha,alpha)));
    //return Beta(2-alpha, alpha)/Beta(2.0-alpha,alpha);
}


double ETw( double alphaA, double alphaAB, double tau ){
    double lambdaA = lambda( alphaA );
    double lambdaAB = lambda( alphaAB );

    return ( 1 - exp( -lambdaA * tau ) ) / lambdaA + exp( -lambdaA * tau) * ( tau + 1 / lambdaAB);
}

double ETb( double alphaAB, double tau ){
    double lambdaAB = lambda( alphaAB );
    return tau + 1/ lambdaAB;
}

double FST_indirect( double alphaA, double alphaAB, double tau ){
    double lambdaA = lambda( alphaA );
    double lambdaAB = lambda( alphaAB );
    return ( 1 - ETw( lambdaA, lambdaAB, tau) / ETb( lambdaAB, tau) );
}

double FST( double alphaA, double alphaAB, double tau ){
    double lambdaA = lambda( alphaA );
    double lambdaAB = lambda( alphaAB );
    cout << "lambdaA = " << lambdaA <<endl;
    return ( 1 - exp( -tau ) ) * ( tau / ( 1 + tau ) );
    //return ( 1 - exp(-lambdaA * tau) ) * ( 1 - 1 / ( tau + 1 / lambdaAB ) / lambdaA );
}

double compute_fst ( vector < valarray <int> > &sites , vector < int > &sample_size ){
    //assert ( sample_size.size() == 2 );
    // at least one mutation in the tree
    // computing within differences
    double Hw = 0, Hb = 0;
    // first add to Hw for population A
    size_t pop_shift = 0;

    for ( size_t pop = 0; pop < sample_size.size(); pop++){
        size_t sample_size_tmp = sample_size[pop];
        double Hw_tmp = 0;
        double Npairs = 0;
        for( size_t sample_i = 0 ; sample_i < (sample_size_tmp-1); sample_i++ ){
            for( size_t sample_j = (sample_i+1); sample_j < sample_size_tmp ; sample_j++ ){
                Npairs++;
                for ( size_t base = 0; base < sites.size(); base++ ){
                    Hw_tmp += abs( sites[base][sample_i+pop_shift] - sites[base][sample_j+pop_shift] );
                }
            }
        }
        //cout << "Differnce in population " <<pop<<" "<< Hw_tmp<<endl;
        Hw += Hw_tmp / ( sample_size_tmp * ( sample_size_tmp - 1) ) ;
        pop_shift += sample_size_tmp;
    }


    double Npairs = 0;
    size_t pop_i_shift = 0;
    size_t pop_j_shift = sample_size[0];
    for ( size_t pop_i = 0; pop_i < (sample_size.size() - 1); pop_i++ ){

        size_t sample_size_pop_i = sample_size[pop_i];

        for( size_t sample_i = 0 ; sample_i < sample_size_pop_i; sample_i++ ){

            for ( size_t pop_j = (pop_i+1); pop_j < sample_size.size() ; pop_j++ ){
                size_t sample_size_pop_j = sample_size[pop_j];
                for( size_t sample_j = 0 ; sample_j < sample_size_pop_j; sample_j++ ){
                    Npairs++;
                    for ( size_t base = 0; base < sites.size(); base++ ){
                        Hb += abs( sites[base][sample_i+pop_i_shift] - sites[base][sample_j+pop_j_shift] ) ;
                    }
                }

            }
        }
        pop_i_shift += sample_size_pop_i;
        pop_j_shift += sample_size[pop_i + 1];
    }
    double size_prod = 1;
    for ( size_t pop = 0; pop < sample_size.size(); pop++ ){
		size_prod *= sample_size[pop];
		}
    Hb /= size_prod;
    //cout << "Hw = " << Hw <<endl;
    //cout << "Hb = " << Hb <<endl;
    //cout << "fst = " << 1.0 - (Hw/Hb) <<endl;

    //return ( 1 -  (Hw*n/(Hb*(n-1)))  );
    return ( 1 -  (Hw / Hb) );
}
