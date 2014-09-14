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

/*! \file sim_gt.cpp
 * \brief Core functions used to simulate gene trees for given species network or tree */

#include "sim_gt.hpp"


action_board::action_board(){
	this->sim_mut_unit_bool  = false;
	this->sim_num_gener_bool = false;
	this->sim_num_mut_bool   = false;
	this->Si_num_bool        = false;
	this->mono_bool          = false;
}

SimulationParameters::SimulationParameters(){
	this->mutation_rate=0.00005;
	//pop_size=10000;
	//pop_size_string_bool=false;
	//mm=2.0;
    //this->net_str = "";
	this->pop_bool=false;
	this->mm_bool=false;
	this->samples_bool=false;
    this->is_Net = false;
	this->num_gener_bool=false;
	this->sp_coal_unit_bool=false;
    this->total_num_lineage = 0;
}


void SimulationParameters::finalize(){
    if ( !this->mm_bool ) { // If coalescent parameter is ungiven, use Kingman coalescent as default
        para_string = write_para_into_tree( net_str, 2.0 ); 
        clog << "Default Kingman coalescent on all branches." << endl;
    }

	if ( !this->pop_bool ) { // If the population size is ungiven, use default population size 10000
        sp_string_pop_size = write_para_into_tree( net_str, 10000.0 ); 
        clog << "Default population size of 10000 on all branches. "<<endl;
    }

	Net net_dummy(net_str);
	
    if ( !net_dummy.is_ultrametric){
	    clog<<"WARNING! NOT ULTRAMETRIC!!!"<<endl;
		//return EXIT_FAIL;
    }
    
    this->is_Net = net_dummy.is_Net;
    
	if ( !this->samples_bool ){
		for (int i = 0; i < net_dummy.tax_name.size(); i++){
			sample_size.push_back(1);
		}
	}
	else{ //  check the number of lineages and the number of species 
		if (sample_size.size() != net_dummy.tax_name.size()) 	throw std::invalid_argument("Numbers of samples and numbers of species not equal!!!");
	}
    
    for ( size_t i = 0; i < this->sample_size.size(); i++) this->total_num_lineage += this->sample_size[i];
    
	sp_string_pop_size = rewrite_pop_string_by_para_string(para_string,sp_string_pop_size);  // checking if modify pop_size_string is needed,

	if ( this->num_gener_bool ) net_str = write_sp_string_in_coal_unit(net_str, sp_string_pop_size);	// Convert number of generations and population size to coalescent unit
	
    sp_string_coal_unit = net_str;
}


sim_one_gt::sim_one_gt ( SimulationParameters* sim_param, action_board* simulation_jobs, ofstream &Si_table): parameters_(sim_param), simulation_jobs_ (simulation_jobs) {
    this->init();
	dout<<"	Starting simulating gene tree from "<<  this->parameters_->sp_string_coal_unit<<endl;
    this->Si_table_ = &Si_table;
	string sp_string_coal_unit= this->parameters_->sp_string_coal_unit;
	string sp_string_pop_size= this->parameters_->sp_string_pop_size;
	string para_string = this->parameters_->para_string;
    
	sim_num_gener_bool_ = this->simulation_jobs_->sim_num_gener_bool;	
    if ( this->simulation_jobs_->sim_mut_unit_bool ||  this->simulation_jobs_->sim_num_mut_bool) sim_num_gener_bool_ = true;
	
    // \todo the following three trees should have been passed with parameters_
    Net my_Net(sp_string_coal_unit);
	Net my_pop_net(sp_string_pop_size);
	Net my_para_net(para_string);
    
	this->initialize_remaining_sp_node ( my_Net );

    this->initialize_gt_tip_nodes( my_Net );    

    this->initialize_gt_internal_nodes ( my_Net.tax_name.size() );

	this->push_in_descdent(); // for removal
     
	int rank_i = 1;
    size_t remaining_sp_node_i=0;	
    
	while ( remaining_sp_node.size() > 0 ){
		size_t node_i = remaining_sp_node[remaining_sp_node_i];
        current_sp_pop_node = &my_Net.NodeContainer[node_i];
		current_sp_pop_size = &my_pop_net.NodeContainer[node_i];
        current_sp_multiCoal_para = &my_para_net.NodeContainer[node_i];
        if ( current_sp_pop_node->rank() == rank_i ){

            this->include_lineages_at_sp_node ( current_sp_pop_node );

			dout << "In population: " << current_sp_pop_node->label << " with rank " << rank_i;
			dout << ((rank_i == my_Net.max_rank) ? ", at the root, everything coalesces":"") ;
            dout << endl;            
			//here to choose go left or right for hybrid node.			
			if ( current_sp_pop_node->hybrid() ) this->assign_lineages_at_sp_node ( current_sp_pop_node );            
			
            double remaining_length = ( rank_i == my_Net.max_rank ) ? 1.0/0.0 : current_sp_pop_node->brchlen1();
            double multi_merge_para = current_sp_multiCoal_para->brchlen1();
            this->implement_coalsecent( current_sp_pop_node->Net_node_contains_gt_node1, remaining_length, multi_merge_para );
            
            if ( rank_i == my_Net.max_rank ) {assert(my_gt_coal_unit.print_all_node_dout());break; }
    
            this->adjust_bl_core (current_sp_pop_node->Net_node_contains_gt_node1, current_sp_pop_node->parent1->height);				
			
			if ( current_sp_pop_node->hybrid() ) {
				dout<<"hybrid node parent1 finished, in parent 2 now"<<endl;
                double remaining_length = current_sp_pop_node->brchlen2();
                double multi_merge_para = current_sp_multiCoal_para->brchlen2();
                this->implement_coalsecent( current_sp_pop_node->Net_node_contains_gt_node2, remaining_length, multi_merge_para);
                this->adjust_bl_core (current_sp_pop_node->Net_node_contains_gt_node2, current_sp_pop_node->parent2->height);
            }
			remaining_sp_node.erase( remaining_sp_node.begin() + remaining_sp_node_i );
		}
		else { 	remaining_sp_node_i++; 	}
	}

    this->finalize( my_Net.tax_name.size() );
}

void sim_one_gt::remove_unused_nodes(){
    for ( size_t node_i = 0; node_i < my_gt_coal_unit.NodeContainer.size(); ){
		if ( my_gt_coal_unit.NodeContainer[node_i].num_descndnt == 0 ){
			my_gt_coal_unit.NodeContainer.erase(my_gt_coal_unit.NodeContainer.begin()+node_i);
			if (sim_num_gener_bool_) my_gt_num_gener.NodeContainer.erase(my_gt_num_gener.NodeContainer.begin() + node_i);
		}
		else{
			node_i++;
		}
	}
}

void sim_one_gt::implement_coalsecent( vector <size_t> & current_alive_lineages, double remaining_length, double multi_merge_para){
    size_t num_lineage = current_alive_lineages.size();
    this->build_lambda_bk_mat( multi_merge_para , num_lineage );			                        
    
    dout << num_lineage << " lineages entering population "<< current_sp_pop_node->label <<" with remaining branch length of " << remaining_length  <<endl;

    // The first condition is valid when proposed coalsecent will happen befoer moving on to the next population. 
    // The second condition is valid when it is above the root
    while ( ( remaining_length > 0) && ( num_lineage > 1) ) {
        size_t new_lineage = remaining_gt_node[0]; // index of the internal node that lineage coalseced to
        this->compute_bl_extension( multi_merge_para, num_lineage );
        dout << "current_lineage_Extension = "<<current_lineage_Extension<< " and remaining_length = "<< remaining_length<<endl;
        if (  current_lineage_Extension < remaining_length ){
            dout << "  " << current_N_lineage_To_Coalesce << " lineages coalesce at time " << current_lineage_Extension << endl;
            for ( size_t nc_i = 0; nc_i < current_N_lineage_To_Coalesce; nc_i++ ){
                size_t lineage_index = rand() % current_alive_lineages.size();
                size_t gt_child_node_index = current_alive_lineages[lineage_index];
                my_gt_coal_unit.NodeContainer[gt_child_node_index].set_brchlen1 ( my_gt_coal_unit.NodeContainer[gt_child_node_index].brchlen1() + current_lineage_Extension ); // b_alpha <- b_alpha + l
                my_gt_coal_unit.NodeContainer[new_lineage].add_child( &my_gt_coal_unit.NodeContainer[gt_child_node_index] );
                my_gt_coal_unit.descndnt[new_lineage] += my_gt_coal_unit.descndnt[gt_child_node_index];
                if (sim_num_gener_bool_){
                    my_gt_num_gener.NodeContainer[gt_child_node_index].set_brchlen1 ( my_gt_num_gener.NodeContainer[gt_child_node_index].brchlen1() + ( current_lineage_Extension * current_sp_pop_size->brchlen1()) ) ; // b_alpha <- b_alpha + l*pop_size
                    my_gt_num_gener.NodeContainer[new_lineage].add_child(&my_gt_num_gener.NodeContainer[gt_child_node_index]);
                    my_gt_num_gener.descndnt[new_lineage] += my_gt_num_gener.descndnt[gt_child_node_index];
                }
                current_alive_lineages.erase(current_alive_lineages.begin()+lineage_index); // X'\{alpha}
            }
            
            my_gt_coal_unit.NodeContainer[new_lineage].num_descndnt = my_gt_coal_unit.descndnt[new_lineage].sum();
            if (sim_num_gener_bool_){
                my_gt_num_gener.NodeContainer[new_lineage].num_descndnt = my_gt_num_gener.descndnt[new_lineage].sum();
            }
            current_alive_lineages.push_back(new_lineage); // introducing a new lineage gamma
            my_gt_coal_unit.NodeContainer[new_lineage].height = my_gt_coal_unit.NodeContainer[new_lineage].child[0]->height + my_gt_coal_unit.NodeContainer[remaining_gt_node[0]].child[0]->brchlen1(); // a_gamma <- a_alpha + b_alpha
            for ( size_t update_brch_i=0;update_brch_i < current_alive_lineages.size()-1;update_brch_i++){
                my_gt_coal_unit.NodeContainer[current_alive_lineages[update_brch_i]].set_brchlen1 ( my_gt_coal_unit.NodeContainer[current_alive_lineages[update_brch_i]].brchlen1()+current_lineage_Extension );		
                if (sim_num_gener_bool_){
                    my_gt_num_gener.NodeContainer[current_alive_lineages[update_brch_i]].set_brchlen1 ( current_lineage_Extension*current_sp_pop_size->brchlen1()+my_gt_num_gener.NodeContainer[current_alive_lineages[update_brch_i]].brchlen1()) ;
                }
            }
            
            remaining_gt_node.erase( remaining_gt_node.begin() );				
            remaining_length -= current_lineage_Extension;
            num_lineage = current_alive_lineages.size();
            dout<<num_lineage<<" live lineages in population "<< current_sp_pop_node->label <<" with remaining branch length of " << remaining_length  <<endl;			
        }
        else return;
    }
}			


void sim_one_gt::finalize( size_t num_taxa ){
    this->remove_unused_nodes();
	//if (!my_gt_coal_unit.check_isUltrametric()){throw "not ultrametic";}	
    this->finalize_gt_str (gt_string_coal_unit , my_gt_coal_unit);
	if (sim_num_gener_bool_) this->finalize_gt_str (gt_string_gener_num , my_gt_num_gener);
        
	// check if the gene tree is ultramatric.
	dout<<"check of if "<<gt_string_coal_unit <<" is ultrametric"<<endl;
	//Net checking_ultra_net(gt_string_coal_unit);
	//if (!checking_ultra_net.is_ultrametric){throw "Gene tree is not ultrametric";}	
	
	if ( this->simulation_jobs_->sim_mut_unit_bool ) build_gt_string_mut_unit( );
	
	if ( this->simulation_jobs_->sim_num_mut_bool ||  this->simulation_jobs_->Si_num_bool ) this->build_mt_tree( );

	if ( num_taxa == 2 && this->simulation_jobs_->mono_bool) compute_monophyly_vec( my_gt_coal_unit, this->parameters_->sample_size );

	dout<<"	end of sim_one_gt::sim_one_gt(sim::param sim_param,action_board my_action)"<<endl;
}


void sim_one_gt::finalize_gt_str( string & gt_tr, Net & gt ){
    gt.NodeContainer.back().CalculateRank();
    gt.rewrite_node_content();
    string gt_tmp_str = gt.NodeContainer.back().node_content + gt.NodeContainer.back().label+";";
    gt_tr = remove_interior_label(gt_tmp_str);
}


void sim_one_gt::build_mt_tree(){
  	double mutation_rate = this->parameters_->mutation_rate;
    Net mt_tree(gt_string_gener_num);
    vector <double> brch_total;
    total_brchlen = 0;
    for ( size_t i = 0 ; i < mt_tree.NodeContainer.size() ; i++ ){
        total_brchlen += mt_tree.NodeContainer[i].brchlen1();
        mt_tree.NodeContainer[i].set_brchlen1( 0.0 );
        brch_total.push_back(total_brchlen);
    }	
    
    brch_total.back() = 0;
    this->total_mut = poisson_rand_var(mutation_rate*total_brchlen);
    for ( int mut_i=0; mut_i < this->total_mut; mut_i++){
        size_t brch_index = 0;
        double u = unifRand()*total_brchlen;
        while ( u > brch_total[brch_index] ){
             brch_index++;
        }
        mt_tree.NodeContainer[brch_index].set_brchlen1 ( mt_tree.NodeContainer[brch_index].brchlen1() + 1 );
    }
    mt_tree.rewrite_node_content();
    gt_string_mut_num = mt_tree.NodeContainer.back().node_content + mt_tree.NodeContainer.back().label + ";";
    gt_string_mut_num = remove_interior_label(gt_string_mut_num);
    
    if ( this->simulation_jobs_->Si_num_bool) this->Si_num_out_table(mt_tree);
}


void sim_one_gt::Si_num_out_table( Net &mt_tree ){
	vector <int> S_i(mt_tree.tip_name.size()-1,0);
	for ( size_t sii = 0; sii < S_i.size(); sii++){
		for ( size_t i = 0; i < mt_tree.NodeContainer.size(); i++){
			int sii_pluse_1 = sii+1;
			if ((mt_tree.NodeContainer[i].brchlen1() > 0) && (sii_pluse_1 == mt_tree.descndnt2[i].sum() )){
				//if ((sii+1)==mt_tree.descndnt2[i].sum()){
					S_i[sii] += mt_tree.NodeContainer[i].brchlen1();
				//}
			}
		}
	}
		
	*Si_table_ << setw(12)<< total_brchlen<<setw(14)<<mt_tree.NodeContainer.back().height<<setw(7)<<this->total_mut<<"  ";
	for ( size_t sii=0;sii<S_i.size();sii++){
		*Si_table_<<setw(4)<<S_i[sii]<<" ";
	}
    *Si_table_ << endl;	
}


void sim_one_gt::adjust_bl_core( vector <size_t> &Net_node_contains_gt_node, double top_time_in_coal_unit){
    dout<<"*************************before adjusting***************"<<endl;
    assert(my_gt_coal_unit.print_all_node_dout());
    for ( size_t i = 0; i < Net_node_contains_gt_node.size(); i++){
        double height_diff_in_coal_unit = top_time_in_coal_unit - my_gt_coal_unit.NodeContainer[Net_node_contains_gt_node[i]].height ;
        if ( sim_num_gener_bool_ ){
            double tmp_bl = ( my_gt_coal_unit.NodeContainer[Net_node_contains_gt_node[i]].height  >  current_sp_pop_node->height ) ? 
                            (height_diff_in_coal_unit) * current_sp_pop_size->brchlen2() : 
                            (height_diff_in_coal_unit - my_gt_coal_unit.NodeContainer[Net_node_contains_gt_node[i]].brchlen1() ) * current_sp_pop_size->brchlen2() + my_gt_num_gener.NodeContainer[Net_node_contains_gt_node[i]].brchlen1() ;
            assert ( (height_diff_in_coal_unit) * current_sp_pop_size->brchlen2() == (height_diff_in_coal_unit - my_gt_coal_unit.NodeContainer[Net_node_contains_gt_node[i]].brchlen1() ) * current_sp_pop_size->brchlen2() + my_gt_num_gener.NodeContainer[Net_node_contains_gt_node[i]].brchlen1());
            //if ( (my_gt_coal_unit.NodeContainer[Net_node_contains_gt_node[i]].height ) >  current_sp_pop_node->height ){
              //my_gt_num_gener.NodeContainer[Net_node_contains_gt_node[i]].set_brchlen1 ( (current_sp_pop_node->parent2->height - my_gt_coal_unit.NodeContainer[Net_node_contains_gt_node[i]].height)*current_sp_pop_size->brchlen2() );//+my_gt_num_gener.NodeContainer[current_sp_pop_node->Net_node_contains_gt_node1[i]]->brchlen1;		
            //}
            //else{
              //my_gt_num_gener.NodeContainer[Net_node_contains_gt_node[i]].set_brchlen1 ( (current_sp_pop_node->parent2->height - my_gt_coal_unit.NodeContainer[Net_node_contains_gt_node[i]].height- my_gt_coal_unit.NodeContainer[Net_node_contains_gt_node[i]].brchlen1() )*current_sp_pop_size->brchlen2() + my_gt_num_gener.NodeContainer[Net_node_contains_gt_node[i]].brchlen1() );
            //}
        }
        my_gt_coal_unit.NodeContainer[Net_node_contains_gt_node[i]].set_brchlen1 ( height_diff_in_coal_unit );
    }
    dout<<"************************* after adjusting***************"<<endl;
    assert(my_gt_coal_unit.print_all_node_dout());
}


void sim_one_gt::build_lambda_bk_mat( double para, size_t num_lineage ){
    lambda_bk_mat.clear();
    if ( para == 2.0 ) return;
    //cout << "para = " << para<<endl;
    assert( para <= 2 );
    assert( para > 0 );
    
	for (size_t b_i = 2;b_i <= num_lineage; b_i++){
		vector <double> lambda_bk_mat_b;
		for (size_t k_i = 2; k_i <= b_i; k_i++){
            // replace this loop by the following
            //double lambda_bk_mat_b_k = ( para < 1 ) ? 
                                    //exp( log( boost::math::binomial_coefficient<double>( b_i, k_i) ) + log( pow( para, (double)k_i ) ) + log( pow ( 1 - para, (double)b_i - (double)k_i ) ) ) :
                                    //exp( log(boost::math::binomial_coefficient<double>( b_i, k_i ) ) + log( Beta( (double)k_i - para, (double)b_i - (double)k_i+para ) ) - log( Beta( 2.0 - para, para ) ) );
            
			double lambda_bk_mat_b_k;
			if ( para < 1 ){//0<psi<1
				//lambda_bk_mat_b_k=n_choose_k(b_i,k_i)*pow(para,k_i)*pow(1-para,b_i-k_i);//.2 is psi lambda_bk=\binom{b}{k}\psi^k (1-\psi)^{b-k}
                //lambda_bk_mat_b_k=exp(log(n_choose_k(b_i,k_i)) + log(pow(para,k_i)) + log(pow(1-para,b_i-k_i)));//.2 is psi lambda_bk=\binom{b}{k}\psi^k (1-\psi)^{b-k}
                //cout<<"normal calculation : "<<exp(log(n_choose_k(b_i,k_i)) + log(pow(para,k_i)) + log(pow(1-para,b_i-k_i)))<<endl;
                lambda_bk_mat_b_k = exp( log( boost::math::binomial_coefficient<double>( b_i, k_i) ) + log( pow( para, (double)k_i ) ) + log( pow ( 1 - para, (double)b_i - (double)k_i ) ) );//.2 is psi lambda_bk=\binom{b}{k}\psi^k (1-\psi)^{b-k}
                //cout<<"boost calculation : "<<lambda_bk_mat_b_k<<endl;
// DEBUG
				//if (isnan(lambda_bk_mat_b_k)){
					//throw std::domain_error("Function \"build_lambda_bk_mat\" returns NaN");
					////dout<<"log(n_choose_k(b_i,k_i)) " <<log(n_choose_k(b_i,k_i))<<endl;
					////dout<<"b="<<b_i<<"  k="<<k_i<<endl;
					////dout<<"log(boost::math::binomial_coefficient(b_i,k_i)) " <<log(boost::math::binomial_coefficient<double>(unsigned(b_i),unsigned(k_i)))<<endl;
					////dout<<"log(pow(para,k_i)) "<< log(pow(para,k_i))<<endl;
					////dout<<"log(pow(1-para,b_i-k_i)) "<<log(pow(1-para,b_i-k_i))<<endl;
				//}
			}
			else{//1<alpha<2
				//lambda_bk_mat_b_k=n_choose_k(b_i,k_i)*Beta(k_i-para,b_i-k_i+para)/Beta(2.0-para,para);// \lambda_{bk}=\binom{b}{k}\frac{B(k-\alpha,b-k+\alpha)}{B(2-\alpha,\alpha)}
                //cout<<n_choose_k(b_i,k_i)*Beta(k_i-para,b_i-k_i+para)/Beta(2.0-para,para)<<endl;
                //lambda_bk_mat_b_k=exp(log(n_choose_k(b_i,k_i))+log(Beta(k_i-para,b_i-k_i+para)) - log(Beta(2.0-para,para)));// \lambda_{bk}=\binom{b}{k}\frac{B(k-\alpha,b-k+\alpha)}{B(2-\alpha,\alpha)}
                //cout<<"normal calculation "<<exp(log(n_choose_k(b_i,k_i))+log(Beta(k_i-para,b_i-k_i+para)) - log(Beta(2.0-para,para)))<<endl;// \lambda_{bk}=\binom{b}{k}\frac{B(k-\alpha,b-k+\alpha)}{B(2-\alpha,\alpha)}
                lambda_bk_mat_b_k = exp( log(boost::math::binomial_coefficient<double>( b_i, k_i ) ) + log( Beta( (double)k_i - para, (double)b_i - (double)k_i+para ) ) - log( Beta( 2.0 - para, para ) ) );// \lambda_{bk}=\binom{b}{k}\frac{B(k-\alpha,b-k+\alpha)}{B(2-\alpha,\alpha)}
				//cout<<"boost result "<<lambda_bk_mat_b_k<<endl;
			}
			lambda_bk_mat_b.push_back(lambda_bk_mat_b_k);
			}
		lambda_bk_mat.push_back(lambda_bk_mat_b);
	}
}


void sim_one_gt::build_gt_string_mut_unit(){
   	double mutation_rate = this->parameters_->mutation_rate;

	Net gt_mut_unit(gt_string_gener_num);
	for ( size_t i = 0; i < gt_mut_unit.NodeContainer.size(); i++){
		gt_mut_unit.NodeContainer[i].set_brchlen1 ( gt_mut_unit.NodeContainer[i].brchlen1() * mutation_rate );
	}
    gt_mut_unit.rewrite_node_content();
    gt_string_mut_unit = gt_mut_unit.NodeContainer.back().node_content + gt_mut_unit.NodeContainer.back().label + ";";
	gt_string_mut_unit=remove_interior_label(gt_string_mut_unit);
}


void sim_one_gt::compute_monophyly_vec( Net &my_gt_coal_unit,vector < int > sample_size){
	vector <double> monophyly_initial(6,0);
	monophyly = monophyly_initial;
	for ( size_t tax_i=0;tax_i<2;tax_i++){
		for ( size_t node_i = 0; node_i < my_gt_coal_unit.NodeContainer.size(); node_i++){
			if (my_gt_coal_unit.descndnt[node_i].sum() == 0){
				break;
			}
			if ((my_gt_coal_unit.descndnt[node_i][tax_i]==sample_size[tax_i]) &&  ( my_gt_coal_unit.descndnt[node_i][tax_i]==my_gt_coal_unit.descndnt[node_i].sum() ) ){
				monophyly[tax_i]=1;
				break;
			}
		}
	}
	if ( monophyly[0] == 1 ){
		if ( monophyly[1] == 1 ) monophyly[2] = 1;
		else                     monophyly[3] = 1;
	}
	else{
		if ( monophyly[1] == 1 ) monophyly[4] = 1;
		else                     monophyly[5] = 1;
	}

	for (size_t mono_i = 0; mono_i < 6; mono_i ++){
		dout<<monophyly[mono_i]<<"  ";
	}
	dout<<endl;
}


/*! \brief Convert the network branch lengths from number of generations to coalescent unit, by dividing the population size 
 * \return string Network in extended newick form, branch lengths are in coalescent unit*/
string SimulationParameters::write_sp_string_in_coal_unit( string &sp_num_gener_string /*! Network in extended newick form, branch lengths are the number of generations*/,
                                                           string &pop_size_string /*! Network in extended newick form, branch lengths are the population sizes*/){
	Net sp_num_gener_net(sp_num_gener_string);
	Net pop_size_net(pop_size_string);
	for ( size_t node_i=0; node_i < sp_num_gener_net.NodeContainer.size(); node_i++){
		if ( node_i < sp_num_gener_net.NodeContainer.size()-1 ){
			sp_num_gener_net.NodeContainer[node_i].set_brchlen1 ( sp_num_gener_net.NodeContainer[node_i].brchlen1() / pop_size_net.NodeContainer[node_i].brchlen1() );
			if ( pop_size_net.NodeContainer[node_i].hybrid() ){
				sp_num_gener_net.NodeContainer[node_i].set_brchlen2 ( sp_num_gener_net.NodeContainer[node_i].brchlen2() / pop_size_net.NodeContainer[node_i].brchlen2() );
			}
		}
	}
    sp_num_gener_net.rewrite_node_content();
	//rewrite_node_content(sp_num_gener_net_node_ptr);
	string sp_coal_unit_string = construct_adding_new_Net_str(sp_num_gener_net);
	
	return sp_coal_unit_string;
}


/*! \brief For alpha coalescent, modify the population size according to the multi merger parameter, N=N^(alpha-1), as mutation must scale with the coalescent timescale, in this case is N^(alpha-1)
 * \return string  Network in extended newick form, branch lengths are the population sizes */
string SimulationParameters::rewrite_pop_string_by_para_string( string para_string /*! Network in extended newick form, branch lengths are the coalescent parameters */,
                                                                string pop_size_string /*! Network in extended newick form, branch lengths are the population sizes*/ ){
	Net para_net_check(para_string);
	Net pop_size_check(pop_size_string);
	for ( size_t node_i = 0; node_i < para_net_check.NodeContainer.size(); node_i++){
		if ( para_net_check.NodeContainer[node_i].brchlen1() < 2 && para_net_check.NodeContainer[node_i].brchlen1() > 1) { // rescale the number of generations for alpha
			pop_size_check.NodeContainer[node_i].set_brchlen1 ( pow( pop_size_check.NodeContainer[node_i].brchlen1(), para_net_check.NodeContainer[node_i].brchlen1() - 1 ) );
		}
        if ( !pop_size_check.NodeContainer[node_i].hybrid() ) continue;
        if ( para_net_check.NodeContainer[node_i].brchlen2() < 2 && para_net_check.NodeContainer[node_i].brchlen2() > 1){
            pop_size_check.NodeContainer[node_i].set_brchlen2 ( pow(pop_size_check.NodeContainer[node_i].brchlen2(), para_net_check.NodeContainer[node_i].brchlen2()-1) );
        }
	}
    pop_size_check.rewrite_node_content();
	string pop_size_string_return = construct_adding_new_Net_str(pop_size_check);
	return pop_size_string_return;
}


double sim_one_gt::update_coal_para( vector < vector <double> > &lambda_bk_mat, double num_lineage ){
	double coal_para = 0;
	for ( int lambda_bk_i = 0; lambda_bk_i <= num_lineage-2; lambda_bk_i++ ){
		coal_para += lambda_bk_mat[num_lineage-2][lambda_bk_i];
		dout<<"lambda_bk_mat[num_lineage-2].size() "<< lambda_bk_mat[num_lineage-2].size() <<endl;
		dout<< " !!! " << num_lineage<<"  "<< lambda_bk_i+2 << endl;
		dout<<"coal_para = "<< coal_para<<endl;
	}
	return coal_para;
}


void sim_one_gt::build_nc_X( size_t num_lineage ){
    this->nc_X = valarray <double> ( num_lineage - 1 );
	for ( size_t kmerge = 0; kmerge < this->nc_X.size(); kmerge++ ){
		this->nc_X[kmerge]= -log( 1-unifRand() ) / lambda_bk_mat[num_lineage-2][kmerge];
		//if (isnan(nc_X[kmerge])){cout<<lambda_bk_mat[num_lineage-2][kmerge]<<endl;}
	}
}

//use heap structure for this!
size_t sim_one_gt::update_nc (){	
	for ( size_t kmerge = 0; kmerge < int(nc_X.size()); kmerge++ ){
		if ( this->nc_X[kmerge] == this->nc_X.min()){
			//nc=kmerge+2;
			return kmerge + 2;
			//break;
		}
	}
    dout << "k merger was never found ... " << endl;
}


void sim_one_gt::init(){
    sim_num_gener_bool_ = false;
}


void sim_one_gt::initialize_gt_tip_nodes( Net & my_Net ){
    dout << " initialize_gt_tip_nodes " << endl;
	for ( size_t i = 0; i < my_Net.NodeContainer.size(); i++ ){
		if ( !my_Net.NodeContainer[i].tip_bool ) continue;
        
        for ( size_t sample_size_i = 0; sample_size_i < this->parameters_->sample_size.size(); sample_size_i++){
            if ( my_Net.descndnt[i][sample_size_i] == 0 ) continue;

            for ( size_t sample_i = 0; sample_i < this->parameters_->sample_size[sample_size_i]; sample_i++){
                this->my_gt_coal_unit.NodeContainer.push_back(my_Net.NodeContainer[i]);
                this->my_gt_coal_unit.descndnt.push_back(my_Net.descndnt[i]);
                this->my_gt_coal_unit.NodeContainer.back().label += "_" + to_string( sample_i+1 );
                this->my_gt_coal_unit.NodeContainer.back().set_brchlen1( 0.0 );
                this->my_gt_coal_unit.NodeContainer.back().parent1 = NULL;
                this->my_gt_coal_unit.NodeContainer.back().parent2 = NULL;
                if ( sim_num_gener_bool_ ){
                    my_gt_num_gener.NodeContainer.push_back(my_gt_coal_unit.NodeContainer.back());
                    my_gt_num_gener.descndnt.push_back(my_gt_coal_unit.descndnt.back());
                }
                assert( my_gt_coal_unit.NodeContainer.back().print_dout( my_Net.is_Net) );
                dout << endl;
                // use the following line to replace mapping_gt_tip_to_sp
                my_Net.NodeContainer[i].Net_node_contains_gt_node1.push_back( this->my_gt_coal_unit.NodeContainer.size()-1 );
            }
        }
        dout << "Population " << my_Net.NodeContainer[i].label<<" have "<<my_Net.NodeContainer[i].Net_node_contains_gt_node1.size() <<" lineages"<<endl;
    }    
}


void sim_one_gt::push_in_descdent(){
    for ( size_t gt_node_i=0; gt_node_i < my_gt_coal_unit.NodeContainer.size(); gt_node_i++){
		valarray <int> descndnt2_dummy;
		my_gt_coal_unit.descndnt2.push_back(descndnt2_dummy);
	}

	if (sim_num_gener_bool_){
		for ( size_t gt_node_i=0;gt_node_i<my_gt_num_gener.NodeContainer.size();gt_node_i++){
			valarray <int> descndnt2_dummy;
			my_gt_num_gener.descndnt2.push_back(descndnt2_dummy);
		}
	}
}


void sim_one_gt::initialize_remaining_sp_node ( Net &my_Net ){
	for ( size_t sp_node_i=0; sp_node_i < my_Net.NodeContainer.size(); sp_node_i++ ) remaining_sp_node.push_back(sp_node_i);
}


void sim_one_gt::initialize_gt_internal_nodes ( size_t num_tax ){
	size_t gt_num_tips = my_gt_coal_unit.NodeContainer.size();
	for ( size_t i = 1; i < gt_num_tips; i++ ){
		Node new_interior_node;
        //new_interior_node.label = ( i == (gt_num_tips-1 ) ) ? "root" : 
                                                                      //"Int_" + to_string(i) ;
		remaining_gt_node.push_back(my_gt_coal_unit.NodeContainer.size()); //my_gt_coal_unit.NodeContainer.size() is the current index of the gt interior node.
		my_gt_coal_unit.NodeContainer.push_back(new_interior_node);
		valarray <int> intialize_descndnt(num_tax);
		my_gt_coal_unit.descndnt.push_back(intialize_descndnt);
		if ( sim_num_gener_bool_ ){
			my_gt_num_gener.NodeContainer.push_back(my_gt_coal_unit.NodeContainer.back());
			my_gt_num_gener.descndnt.push_back(my_gt_coal_unit.descndnt.back());
		}
	}
}
	
void sim_one_gt::include_lineages_at_sp_node( Node * sp_node ){
    for ( size_t i = 0; i < sp_node->child.size(); i++){
        //if (sp_node->child[i]->parent1->label == sp_node->label){
        if ( sp_node->child[i]->parent1 == sp_node ){ // to be used, check later
            for ( size_t j = 0; j < sp_node->child[i]->Net_node_contains_gt_node1.size(); j++){
                sp_node->Net_node_contains_gt_node1.push_back ( sp_node->child[i]->Net_node_contains_gt_node1[j]);
            }
        }
        else if ( sp_node->child[i]->parent2 != sp_node ){
            assert ( sp_node->child[i]->parent2 != sp_node );
            for ( size_t j = 0; j < sp_node->child[i]->Net_node_contains_gt_node2.size(); j++){
                sp_node->Net_node_contains_gt_node1.push_back ( sp_node->child[i]->Net_node_contains_gt_node2[j]);
            }					
        }
        else continue;
    }
}


void sim_one_gt::assign_lineages_at_sp_node ( Node *sp_node ){
    dout<<"hybrid node"<<endl;
    vector < size_t > index_container_tmp = sp_node->Net_node_contains_gt_node1;
    sp_node->Net_node_contains_gt_node1.clear();
    double left_para = sp_node->extract_hybrid_para();
    for ( size_t i = 0; i < index_container_tmp.size(); i++ ){
        if ( unifRand() < left_para ) sp_node->Net_node_contains_gt_node1.push_back(index_container_tmp[i]);
        else                          sp_node->Net_node_contains_gt_node2.push_back(index_container_tmp[i]);
    }
}


void sim_one_gt::compute_bl_extension( double multi_merge_para, size_t num_lineage ){
    current_N_lineage_To_Coalesce = 0; // reset
    current_lineage_Extension = 0; // reset
    if ( num_lineage <= 1) return;

    if ( multi_merge_para != 2 ) this->build_nc_X( num_lineage );
    current_N_lineage_To_Coalesce = (multi_merge_para == 2) ? 2 : update_nc();
    current_lineage_Extension = (multi_merge_para == 2) ? kingman_bl(num_lineage) : current_lineage_Extension = nc_X.min();
    
    assert(current_N_lineage_To_Coalesce > 1);
}


/*! \brief Write a fixed parameter into a externed newick formatted network string*/
string write_para_into_tree( string in_str /*! Externed newick formatted network string*/, 
                             double para /*! Coalescent parameter or fixed population sizes */){
	if ( in_str.size() == 0 ) throw std::invalid_argument("Please define the input tree (network).");

	Net para_Net(in_str);
	for ( size_t i = 0; i < para_Net.NodeContainer.size(); i++){
		para_Net.NodeContainer[i].set_brchlen1( para );
		if ( para_Net.NodeContainer[i].hybrid() )  para_Net.NodeContainer[i].set_brchlen2( para );
	}
    para_Net.rewrite_node_content();
	return construct_adding_new_Net_str(para_Net);
}


string construct_adding_new_Net_str( Net &in_Net ){
    string out_str = in_Net.NodeContainer.back().node_content + in_Net.NodeContainer.back().label;
	if ( in_Net.NodeContainer.back().brchlen1() != 0 ) 	out_str += ":" + to_string( in_Net.NodeContainer.back().brchlen1() );
	out_str.push_back(';');
	return out_str;
}
