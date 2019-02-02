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

/*! \file sim_gt.cpp
 * \brief Core functions used to simulate gene trees for given species network or tree */

#include "sim_gt.hpp"

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
        for (size_t i = 0; i < net_dummy.tax_name.size(); i++){
            sample_size.push_back(1);
        }
    }
    else{ //  check the number of lineages and the number of species
        if (sample_size.size() != net_dummy.tax_name.size())     throw std::invalid_argument("Numbers of samples and numbers of species not equal!!!");
    }

    for ( size_t i = 0; i < this->sample_size.size(); i++) this->total_num_lineage += this->sample_size[i];

    sp_string_pop_size = rewrite_pop_string_by_para_string(para_string, sp_string_pop_size);  // checking if modify pop_size_string is needed,

    if ( this->num_gener_bool ) {
        net_str = write_sp_string_in_coal_unit(net_str, sp_string_pop_size);    // Convert number of generations and population size to coalescent unit
    }

    sp_string_coal_unit = net_str;
    my_Net = new Net (sp_string_coal_unit);
    my_pop_net = new Net (sp_string_pop_size);
    my_para_net = new Net (para_string);
}




//simTree::simTree ( SimulationParameters* sim_param, action_board* simulation_jobs ): parameters_(sim_param), simulation_jobs_ (simulation_jobs) {
    //this->init();
    //dout<<"    Starting simulating gene tree from "<<  this->parameters_->sp_string_coal_unit<<endl;
    //this->core();
//}


simTree::simTree ( SimulationParameters* sim_param, action_board* simulation_jobs, ofstream &Si_table, MersenneTwister *mt ): parameters_(sim_param), simulation_jobs_ (simulation_jobs) {
    this->mt = mt;
    this->init();
    dout<<"    Starting simulating gene tree from "<<  this->parameters_->sp_string_coal_unit<<endl;
    this->Si_table_ = &Si_table;
    this->core();
}

void simTree::core (){
    sim_num_gener_bool_ = this->simulation_jobs_->sim_num_gener_bool;
    if ( this->simulation_jobs_->sim_mut_unit_bool ||  this->simulation_jobs_->sim_num_mut_bool) {
        sim_num_gener_bool_ = true;
    }

    this->initialize_remaining_sp_node ( this->parameters_->my_Net );
    this->initialize_gt_tip_nodes( this->parameters_->my_Net );
    this->initialize_gt_internal_nodes ( this->parameters_->my_Net->tax_name.size() );
    this->push_in_descdent(); // for removal
    int rank_i = 1;
    size_t remaining_sp_node_i=0;

    while ( remaining_sp_node.size() > 0 ){
        size_t node_i = remaining_sp_node[remaining_sp_node_i];
        if ( this->parameters_->my_Net->NodeContainer[node_i].rank() == rank_i ){
            //this->include_lineages_at_sp_node ( current_sp_pop_node );
            this->include_lineages_at_sp_node ( this->parameters_->my_Net->NodeContainer[node_i] );

            dout << "In population: " << this->parameters_->my_Net->NodeContainer[node_i].label << " with rank " << rank_i;
            dout << ((rank_i == this->parameters_->my_Net->max_rank) ? ", at the root, everything coalesces":"") ;
            dout << endl;

            //here to choose go left or right for hybrid node.
            if ( this->parameters_->my_Net->NodeContainer[node_i].hybrid() ) {
                this->assign_lineages_at_sp_node ( this->parameters_->my_Net->NodeContainer[node_i] );
            }

            //
            double remaining_length = ( rank_i == this->parameters_->my_Net->max_rank ) ? 1.0/0.0 : this->parameters_->my_Net->NodeContainer[node_i].brchlen1();
            double multi_merge_para = this->parameters_->my_para_net->NodeContainer[node_i].brchlen1();
            double pop_size = this->parameters_->my_pop_net->NodeContainer[node_i].brchlen1();

            remaining_length = this->implement_coalsecent( this->parameters_->my_Net->NodeContainer[node_i].Net_node_contains_gt_node1,
                                        remaining_length, multi_merge_para, pop_size,
                                        this->parameters_->my_Net->NodeContainer[node_i].label);

            if ( rank_i == this->parameters_->my_Net->max_rank ) {
                assert(this->print_all_node_dout());
                if (sim_num_gener_bool_){
                    assert(my_gt_num_gener.print_all_node_dout());
                }
                break;
            }

            this->adjust_bl_core ( this->parameters_->my_Net->NodeContainer[node_i].Net_node_contains_gt_node1,
                                   this->parameters_->my_Net->NodeContainer[node_i].height(),
                                   this->parameters_->my_Net->NodeContainer[node_i].brchlen1(),
                                   pop_size, remaining_length );

            if ( this->parameters_->my_Net->NodeContainer[node_i].hybrid() ) {
                dout<<"hybrid node parent1 finished, in parent 2 now"<<endl;
                double remaining_length = this->parameters_->my_Net->NodeContainer[node_i].brchlen2();
                double multi_merge_para = this->parameters_->my_para_net->NodeContainer[node_i].brchlen2();
                double pop_size = this->parameters_->my_pop_net->NodeContainer[node_i].brchlen2();
                remaining_length = this->implement_coalsecent( this->parameters_->my_Net->NodeContainer[node_i].Net_node_contains_gt_node2,
                                            remaining_length, multi_merge_para, pop_size,
                                            this->parameters_->my_Net->NodeContainer[node_i].label);
                this->adjust_bl_core ( this->parameters_->my_Net->NodeContainer[node_i].Net_node_contains_gt_node2,
                                       this->parameters_->my_Net->NodeContainer[node_i].height(),
                                       this->parameters_->my_Net->NodeContainer[node_i].brchlen2(),
                                       pop_size, remaining_length);

            }
            remaining_sp_node.erase( remaining_sp_node.begin() + remaining_sp_node_i );
        }
        else {     remaining_sp_node_i++;     }

        if ( remaining_sp_node_i == remaining_sp_node.size() - 1 ){
            rank_i++;
            remaining_sp_node_i=0;
        }
    }

    this->finalize( this->parameters_->my_Net->tax_name.size() );
}


void simTree::remove_unused_nodes(){
    for ( size_t node_i = 0; node_i < this->NodeContainer.size(); ){
        if ( this->NodeContainer[node_i].num_descndnt == 0 ){
            this->NodeContainer.erase(this->NodeContainer.begin()+node_i);
            if (sim_num_gener_bool_) {
                my_gt_num_gener.NodeContainer.erase(my_gt_num_gener.NodeContainer.begin() + node_i);
            }
        }
        else{
            node_i++;
        }
    }
}

double simTree::implement_coalsecent( vector <size_t> & current_alive_lineages,
                                    double remaining_length, double multi_merge_para, double pop_size, string node_label){
    this->current_lineage_Extension = 0;
    size_t num_lineage = current_alive_lineages.size();
    this->build_lambda_bk_mat( multi_merge_para , num_lineage );

    dout << num_lineage << " lineages entering population "<< node_label <<" with remaining branch length of " << remaining_length  <<endl;

    // The first condition is valid when proposed coalsecent will happen befoer moving on to the next population.
    // The second condition is valid when it is above the root
    while ( ( remaining_length > 0) && ( num_lineage > 1) ) {
        size_t new_lineage = remaining_gt_node[0]; // index of the internal node that lineage coalseced to
        this->compute_bl_extension( multi_merge_para, num_lineage );
        dout << "current_lineage_Extension = "<<current_lineage_Extension<< " and remaining_length = "<< remaining_length<<endl;
        if (sim_num_gener_bool_){
            dout << "current_lineage_Extension in generation = "<<current_lineage_Extension*pop_size<< " and remaining_length in generation = "<< remaining_length*pop_size<<endl;
        }
        if (  current_lineage_Extension < remaining_length ){
            dout << "  " << current_N_lineage_To_Coalesce << " lineages coalesce at time " << current_lineage_Extension << endl;
            for ( size_t nc_i = 0; nc_i < current_N_lineage_To_Coalesce; nc_i++ ){

                // DEBUG, problem is here!!!
                //MTRand_int32
                size_t lineage_index = mt->sampleInt(current_alive_lineages.size());
                //size_t lineage_index = rand() % current_alive_lineages.size();
                // DEBUG
                size_t gt_child_node_index = current_alive_lineages[lineage_index];
                this->NodeContainer[gt_child_node_index].set_brchlen1 ( this->NodeContainer[gt_child_node_index].brchlen1() + current_lineage_Extension ); // b_alpha <- b_alpha + l
                this->NodeContainer[new_lineage].add_child( &this->NodeContainer[gt_child_node_index] );
                this->descndnt[new_lineage] += this->descndnt[gt_child_node_index];
                if (sim_num_gener_bool_){
                    my_gt_num_gener.NodeContainer[gt_child_node_index].set_brchlen1 ( my_gt_num_gener.NodeContainer[gt_child_node_index].brchlen1() + ( current_lineage_Extension * pop_size) ) ; // b_alpha <- b_alpha + l*pop_size
                    my_gt_num_gener.NodeContainer[new_lineage].add_child(&my_gt_num_gener.NodeContainer[gt_child_node_index]);
                    my_gt_num_gener.descndnt[new_lineage] += my_gt_num_gener.descndnt[gt_child_node_index];
                }
                current_alive_lineages.erase(current_alive_lineages.begin()+lineage_index); // X'\{alpha}
            }

            this->NodeContainer[new_lineage].num_descndnt = this->descndnt[new_lineage].sum();
            if (sim_num_gener_bool_){
                my_gt_num_gener.NodeContainer[new_lineage].num_descndnt = my_gt_num_gener.descndnt[new_lineage].sum();
            }
            current_alive_lineages.push_back(new_lineage); // introducing a new lineage gamma
            this->NodeContainer[new_lineage].set_height( this->NodeContainer[new_lineage].child[0]->height() + this->NodeContainer[remaining_gt_node[0]].child[0]->brchlen1() ); // a_gamma <- a_alpha + b_alpha
            if (sim_num_gener_bool_){
                my_gt_num_gener.NodeContainer[new_lineage].set_height( my_gt_num_gener.NodeContainer[new_lineage].child[0]->height() + my_gt_num_gener.NodeContainer[remaining_gt_node[0]].child[0]->brchlen1() ); // a_gamma <- a_alpha + b_alpha
                my_gt_num_gener.NodeContainer[new_lineage].set_brchlen1(0); // a_gamma <- a_alpha + b_alpha
            }

            for ( size_t update_brch_i=0;update_brch_i < current_alive_lineages.size()-1;update_brch_i++){
                this->NodeContainer[current_alive_lineages[update_brch_i]].set_brchlen1 ( this->NodeContainer[current_alive_lineages[update_brch_i]].brchlen1()+current_lineage_Extension );
                if (sim_num_gener_bool_){
                    my_gt_num_gener.NodeContainer[current_alive_lineages[update_brch_i]].set_brchlen1 ( current_lineage_Extension*pop_size + my_gt_num_gener.NodeContainer[current_alive_lineages[update_brch_i]].brchlen1()) ;
                }
            }

            remaining_gt_node.erase( remaining_gt_node.begin() );
            remaining_length -= current_lineage_Extension;
            num_lineage = current_alive_lineages.size();
            dout<<num_lineage<<" live lineages in population "<< node_label <<" with remaining branch length of " << remaining_length  <<endl;
        } else {
            dout << " Nothing to coalesce. " << endl;
            current_lineage_Extension = 0;
            return remaining_length;
        }
    }
    return remaining_length;
}


void simTree::finalize( size_t num_taxa ){
    this->remove_unused_nodes();

    this->NodeContainer.back().CalculateRank();
    this->rewrite_node_content();
    string gt_tmp_str = this->NodeContainer.back().node_content + this->NodeContainer.back().label+";";
    gt_string_coal_unit = remove_interior_label(gt_tmp_str);
    //gt_string_coal_unit = print_newick( & this->NodeContainer.back() ) + ";";
    //this->finalize_gt_str (gt_string_coal_unit , (*this) );
    if (sim_num_gener_bool_) {
        this->finalize_gt_str (gt_string_gener_num , my_gt_num_gener);
    }

    // check if the gene tree is ultramatric.
    dout<<"check of if "<<gt_string_coal_unit <<" is ultrametric"<<endl;
    this->check_isUltrametric();
    //if ( !this->is_ultrametric ){ throw std::invalid_argument ("Gene tree is not ultrametric"); }
    if ( !this->is_ultrametric ){ clog << "WARNING: Gene tree is not ultrametric" << endl; }

    if ( this->simulation_jobs_->sim_mut_unit_bool ) { build_gt_string_mut_unit( ); dout<<"mut_tree built"<<endl;}

    if ( this->simulation_jobs_->sim_num_mut_bool ||  this->simulation_jobs_->Si_num_bool ) { this->build_mt_tree( );   dout<<"mt_tree built"<<endl; }

    if ( num_taxa == 2 && this->simulation_jobs_->mono_bool) compute_monophyly_vec( this->parameters_->sample_size );

    dout<<"    end of simTree::simTree(sim::param sim_param,action_board my_action)"<<endl;
}


void simTree::finalize_gt_str( string & gt_str, Tree & gt ){
    //gt.NodeContainer.back().CalculateRank();
    //gt.rewrite_node_content();
    //string gt_tmp_str = gt.NodeContainer.back().node_content + gt.NodeContainer.back().label+";";
    //gt_str = remove_interior_label(gt_tmp_str);
    gt_str = print_newick( &gt.NodeContainer.back() ) + ";";
}


void simTree::build_mt_tree(){
    double mutation_rate = this->parameters_->mutation_rate;
    Tree mt_tree(gt_string_gener_num);
    vector <double> brch_total;
    total_brchlen = 0;
    for ( size_t i = 0 ; i < mt_tree.NodeContainer.size() ; i++ ){
        total_brchlen += mt_tree.NodeContainer[i].brchlen1();
        mt_tree.NodeContainer[i].set_brchlen1( 0.0 );
        brch_total.push_back(total_brchlen);
    }

    brch_total.back() = 0;
    this->total_mut = poisson_rand_var(mutation_rate*total_brchlen);
    for ( size_t mut_i = 0; mut_i < this->total_mut; mut_i++){
        size_t brch_index = 0;
        double u = unifRand()*total_brchlen;
        while ( u > brch_total[brch_index] ){
             brch_index++;
        }
        mt_tree.NodeContainer[brch_index].set_brchlen1 ( mt_tree.NodeContainer[brch_index].brchlen1() + 1 );
    }
    //mt_tree.rewrite_node_content();
    //gt_string_mut_num = mt_tree.NodeContainer.back().node_content + mt_tree.NodeContainer.back().label + ";";
    //gt_string_mut_num = remove_interior_label(gt_string_mut_num);
    gt_string_mut_num = mt_tree.print_newick ( & mt_tree.NodeContainer.back() ) + ";";
    if ( this->simulation_jobs_->Si_num_bool) this->Si_num_out_table(mt_tree);
}


void simTree::Si_num_out_table( Tree &mt_tree ){
    vector <int> S_i(mt_tree.tip_name.size()-1,0);
    for ( size_t sii = 0; sii < S_i.size(); sii++){
        for ( size_t i = 0; i < mt_tree.NodeContainer.size(); i++){
            int sii_pluse_1 = sii+1;
            if ((mt_tree.NodeContainer[i].brchlen1() > 0) && (sii_pluse_1 == mt_tree.samples_below[i].sum() )){
                //if ((sii+1)==mt_tree.samples_below[i].sum()){
                    S_i[sii] += mt_tree.NodeContainer[i].brchlen1();
                //}
            }
        }
    }

    *Si_table_ << total_brchlen<<'\t'<<mt_tree.NodeContainer.back().height() <<'\t'<<this->total_mut;
    for ( size_t sii=0;sii<S_i.size();sii++){
        *Si_table_<<'\t'<<S_i[sii];
    }
    *Si_table_ << endl;
}


void simTree::adjust_bl_core( vector <size_t> &Net_node_contains_gt_node, double bottom_time_in_coal_unit, double pop_bl_in_coal_unit, double pop_size, double remaining_length ){
    double top_time_in_coal_unit = bottom_time_in_coal_unit + pop_bl_in_coal_unit;
    dout<<"*************************before adjusting***************"<<endl;
    assert(this->print_all_node_dout());
    assert(this->my_gt_num_gener.print_all_node_dout());
    for ( size_t i = 0; i < Net_node_contains_gt_node.size(); i++){
        double height_diff_in_coal_unit = top_time_in_coal_unit - this->NodeContainer[Net_node_contains_gt_node[i]].height() ;
        this->NodeContainer[Net_node_contains_gt_node[i]].set_brchlen1(this->NodeContainer[Net_node_contains_gt_node[i]].brchlen1()+remaining_length);
        if ( sim_num_gener_bool_ ){
            double current_node_bl_in_num_gener = my_gt_num_gener.NodeContainer[Net_node_contains_gt_node[i]].brchlen1();
            if (pop_bl_in_coal_unit>current_lineage_Extension & current_lineage_Extension > 0){
                my_gt_num_gener.NodeContainer[Net_node_contains_gt_node[i]].set_brchlen1 ( current_node_bl_in_num_gener + (pop_bl_in_coal_unit-current_lineage_Extension ) * pop_size );
            } else {
                my_gt_num_gener.NodeContainer[Net_node_contains_gt_node[i]].set_brchlen1 ( current_node_bl_in_num_gener + (pop_bl_in_coal_unit) * pop_size );
            }
        }
    }
    dout<<"************************* after adjusting***************"<<endl;
    assert(this->print_all_node_dout());
    if ( sim_num_gener_bool_ ){
        assert(my_gt_num_gener.print_all_node_dout());
    }
}


void simTree::build_lambda_bk_mat( double para, size_t num_lineage ){
    if ( para > 2 or para < 0){ // TODO: move this check to early on?
        throw std::invalid_argument(string("Multiple merger parameter ") + to_string (static_cast<double>(para)) + string(" is out of the range of [0, 2]."));
    }
    assert( para <= 2 );
    assert( para >= 0 );
    lambda_bk_mat.clear();
    if ( para == 2.0 ) return;

    assert( para < 2 );
    for ( int b_int = 2; b_int <= num_lineage; b_int++){
        double b_i = (double)b_int;
        vector <double> lambda_bk_mat_b;
        for ( int k_int = 2; k_int <= b_int; k_int++){
            double lambda_bk_mat_b_k = ( para > 1 ) ? this->lambdaAlpha( b_i, (double)k_int, para):
                                                      this->lambdaPsi( b_i, (double)k_int, para);
            lambda_bk_mat_b.push_back(lambda_bk_mat_b_k);
        }
        lambda_bk_mat.push_back(lambda_bk_mat_b);
    }
}

long double logbinomial(long double n, long double k) {
  assert(n >= 0.0L);
  assert(k >= 0.0L);
  assert(k <= n);

  long double x= 0.0L;
  long double bi= 0.0L;

  while (k > x)
  {
    bi+= logl(n-x);
    bi-= logl(k-x);
    x++;
  }
  return bi;
}

double simTree::lambdaAlpha( double b, double k, double para ){
    assert ( b >= k );
    assert ( k > 1 );
    //if ( b < k) throw std::invalid_argument("b can not be less than k");
    //  \lambda_{bk}=\binom{b}{k}\frac{B(k-\alpha,b-k+\alpha)}{B(2-\alpha,\alpha)}
    return exp( log( binomial_coefficient( b, k)) + log( Beta ( k - para, b-k + para) ) - log( Beta( 2.0 - para, para ) ) );
}

double simTree::lambdaPsi( double b, double k, double para ){
    assert ( b >= k );
    assert ( k > 1.0 );
    //if ( b < k) throw std::invalid_argument("b can not be less than k");
    //.2 is psi lambda_bk=\binom{b}{k}\psi^{k-2} (1-\psi)^{b-k}
    if ( para == 1.0 && (int)b == (int)k ){
        return 1.0;
    } else if ( para == 1.0 && (int)b != (int)k ){
        return 0.0;
    } else if ( para == 0.0 && k == 2.0 ){
        return binomial_coefficient(b, k);
    } else if ( para == 0.0 && k != 2.0 ){
        return 0.0;
    } else {
        return exp ( logbinomial(b,k) + (k-2)*log(para) + (b-k)*log(1.0-para) );
    }
}

void simTree::build_gt_string_mut_unit(){
    double mutation_rate = this->parameters_->mutation_rate;

    Net gt_mut_unit(gt_string_gener_num);
    for ( size_t i = 0; i < gt_mut_unit.NodeContainer.size(); i++){
        gt_mut_unit.NodeContainer[i].set_brchlen1 ( gt_mut_unit.NodeContainer[i].brchlen1() * mutation_rate );
    }
    //gt_mut_unit.rewrite_node_content();
    //gt_string_mut_unit = gt_mut_unit.NodeContainer.back().node_content + gt_mut_unit.NodeContainer.back().label + ";";
    //gt_string_mut_unit=remove_interior_label(gt_string_mut_unit);
    gt_string_mut_unit = gt_mut_unit.print_newick( &gt_mut_unit.NodeContainer.back() ) + ";";
}


void simTree::compute_monophyly_vec( vector < int > sample_size ){
    vector <double> monophyly_initial(6,0);
    monophyly = monophyly_initial;
    for ( size_t tax_i=0;tax_i<2;tax_i++){
        for ( size_t node_i = 0; node_i < this->NodeContainer.size(); node_i++){
            if (this->descndnt[node_i].sum() == 0){
                break;
            }
            if ((this->descndnt[node_i][tax_i]==sample_size[tax_i]) &&  ( this->descndnt[node_i][tax_i]==this->descndnt[node_i].sum() ) ){
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

    bool is_equal = false;
    if(sp_num_gener_net.tax_name.size() == pop_size_net.tax_name.size()){
      is_equal = std::equal(sp_num_gener_net.tax_name.begin(), sp_num_gener_net.tax_name.end(), pop_size_net.tax_name.begin());
    }

    if (!is_equal){
        throw std::invalid_argument("Species tree taxa name do not match!");
    }

    for ( size_t node_i=0; node_i < sp_num_gener_net.NodeContainer.size(); node_i++){
        if ( node_i < sp_num_gener_net.NodeContainer.size()-1 ){
            sp_num_gener_net.NodeContainer[node_i].set_brchlen1 ( sp_num_gener_net.NodeContainer[node_i].brchlen1() / pop_size_net.NodeContainer[node_i].brchlen1() );
            if ( pop_size_net.NodeContainer[node_i].hybrid() ){
                sp_num_gener_net.NodeContainer[node_i].set_brchlen2 ( sp_num_gener_net.NodeContainer[node_i].brchlen2() / pop_size_net.NodeContainer[node_i].brchlen2() );
            }
        }
    }
    sp_num_gener_net.rewrite_node_content();
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


double simTree::update_coal_para( vector < vector <double> > &lambda_bk_mat, double num_lineage ){
    double coal_para = 0;
    for ( int lambda_bk_i = 0; lambda_bk_i <= num_lineage-2; lambda_bk_i++ ){
        coal_para += lambda_bk_mat[num_lineage-2][lambda_bk_i];
        dout<<"lambda_bk_mat[num_lineage-2].size() "<< lambda_bk_mat[num_lineage-2].size() <<endl;
        dout<< " !!! " << num_lineage<<"  "<< lambda_bk_i+2 << endl;
        dout<<"coal_para = "<< coal_para<<endl;
    }
    return coal_para;
}


void simTree::build_nc_X( size_t num_lineage ){
    this->nc_X = valarray <double> ( num_lineage - 1 );
    for ( size_t kmerge = 0; kmerge < this->nc_X.size(); kmerge++ ){
        this->nc_X[kmerge]= -log( 1-unifRand() ) / lambda_bk_mat[num_lineage-2][kmerge];
        //if (isnan(nc_X[kmerge])){cout<<lambda_bk_mat[num_lineage-2][kmerge]<<endl;}
    }
}

//use heap structure for this!
size_t simTree::update_nc (){
    size_t kmerge = 0;
    for ( ; kmerge < int(nc_X.size()); kmerge++ ){
        if ( this->nc_X[kmerge] == this->nc_X.min()){
            //nc=kmerge+2;

            break;
        }
    }
    assert ( this->nc_X[kmerge] == this->nc_X.min() );
    if ( this->nc_X[kmerge] != this->nc_X.min() )
        throw std::invalid_argument ( "k merger was never found ... " );
    return kmerge + 2;
}


void simTree::initialize_gt_tip_nodes( Net * my_Net ){
    dout << " initialize_gt_tip_nodes " << endl;
    for ( size_t i = 0; i < my_Net->NodeContainer.size(); i++ ){
        my_Net->NodeContainer[i].Net_node_contains_gt_node1.clear();
        my_Net->NodeContainer[i].Net_node_contains_gt_node2.clear();
        if ( !my_Net->NodeContainer[i].tip_bool ) continue;

        for ( size_t sample_size_i = 0; sample_size_i < this->parameters_->sample_size.size(); sample_size_i++){
            if ( my_Net->descndnt[i][sample_size_i] == 0 ) continue;

            for ( size_t sample_i = 0; sample_i < this->parameters_->sample_size[sample_size_i]; sample_i++){
                this->NodeContainer.push_back(my_Net->NodeContainer[i]);
                this->descndnt.push_back(my_Net->descndnt[i]);
                this->NodeContainer.back().label += "_" + to_string(static_cast<long long>(sample_i+1 ));
                this->NodeContainer.back().set_brchlen1( 0.0 );
                this->NodeContainer.back().parent1 = NULL;
                this->NodeContainer.back().parent2 = NULL;
                if ( sim_num_gener_bool_ ){
                    my_gt_num_gener.NodeContainer.push_back(this->NodeContainer.back());
                    my_gt_num_gener.descndnt.push_back(this->descndnt.back());
                }
                assert( this->NodeContainer.back().print_dout( my_Net->is_Net) );
                dout << endl;
                // use the following line to replace mapping_gt_tip_to_sp
                my_Net->NodeContainer[i].Net_node_contains_gt_node1.push_back( this->NodeContainer.size()-1 );
            }
        }
        dout << "Population " << my_Net->NodeContainer[i].label<<" have "<<my_Net->NodeContainer[i].Net_node_contains_gt_node1.size() <<" lineages"<<endl;
    }
}


void simTree::push_in_descdent(){
    for ( size_t gt_node_i=0; gt_node_i < this->NodeContainer.size(); gt_node_i++){
        valarray <int> samples_below_dummy;
        this->samples_below.push_back(samples_below_dummy);
    }

    if (sim_num_gener_bool_){
        for ( size_t gt_node_i=0;gt_node_i<my_gt_num_gener.NodeContainer.size();gt_node_i++){
            valarray <int> samples_below_dummy;
            my_gt_num_gener.samples_below.push_back(samples_below_dummy);
        }
    }
}


void simTree::initialize_remaining_sp_node ( Net *my_Net ){
    for ( size_t sp_node_i=0; sp_node_i < my_Net->NodeContainer.size(); sp_node_i++ ) remaining_sp_node.push_back(sp_node_i);
}


void simTree::initialize_gt_internal_nodes ( size_t num_tax ){
    size_t gt_num_tips = this->NodeContainer.size();
    for ( size_t i = 1; i < gt_num_tips; i++ ){
        Node new_interior_node;
        //new_interior_node.label = ( i == (gt_num_tips-1 ) ) ? "root" :
                                                                      //"Int_" + to_string(i) ;
        remaining_gt_node.push_back(this->NodeContainer.size()); //this->NodeContainer.size() is the current index of the gt interior node.
        this->NodeContainer.push_back(new_interior_node);
        valarray <int> intialize_descndnt(num_tax);
        this->descndnt.push_back(intialize_descndnt);
        if ( sim_num_gener_bool_ ){
            my_gt_num_gener.NodeContainer.push_back(this->NodeContainer.back());
            my_gt_num_gener.descndnt.push_back(this->descndnt.back());
        }
    }
}

void simTree::include_lineages_at_sp_node( Node & sp_node ){
    dout << "include lineages at sp_node " << &sp_node << " " << sp_node.label << endl;
    dout << "include lineage: ";
    for ( size_t i = 0; i < sp_node.child.size(); i++){
        if ( sp_node.child[i]->parent1 == &sp_node ){ // to be used, check later
            for ( size_t j = 0; j < sp_node.child[i]->Net_node_contains_gt_node1.size(); j++){
                bool unique = true;
                for (size_t check_i = 0; check_i < sp_node.Net_node_contains_gt_node1.size(); check_i++) {
                    if (sp_node.child[i]->Net_node_contains_gt_node1[j] == sp_node.Net_node_contains_gt_node1[check_i]){
                        unique = false;
                        break;
                    }
                }

                if (unique){
                    sp_node.Net_node_contains_gt_node1.push_back ( sp_node.child[i]->Net_node_contains_gt_node1[j]);
                    dout << &this->NodeContainer[sp_node.Net_node_contains_gt_node1.back()] << ", ";
                }
            }
            if (sp_node.child[i]->parent2 == &sp_node ){ // This is for looping case, 1 sp node has two path to a hybrid node
                for ( size_t j = 0; j < sp_node.child[i]->Net_node_contains_gt_node2.size(); j++){
                    bool unique = true;
                    for (size_t check_i = 0; check_i < sp_node.Net_node_contains_gt_node1.size(); check_i++) {
                        if (sp_node.child[i]->Net_node_contains_gt_node2[j] == sp_node.Net_node_contains_gt_node1[check_i]){
                            unique = false;
                            break;
                        }
                    }

                    if (unique){
                        sp_node.Net_node_contains_gt_node1.push_back ( sp_node.child[i]->Net_node_contains_gt_node2[j]);
                        dout << &this->NodeContainer[sp_node.Net_node_contains_gt_node1.back()] << ", ";
                    }
                }
            }
        } else if ( sp_node.child[i]->parent2 == &sp_node ){
            assert ( sp_node.child[i]->parent2 == &sp_node );
            for ( size_t j = 0; j < sp_node.child[i]->Net_node_contains_gt_node2.size(); j++){
                bool unique = true;
                for (size_t check_i = 0; check_i < sp_node.Net_node_contains_gt_node1.size(); check_i++) {
                    if (sp_node.child[i]->Net_node_contains_gt_node2[j] == sp_node.Net_node_contains_gt_node1[check_i]){
                        unique = false;
                        break;
                    }
                }
                if (unique){
                    sp_node.Net_node_contains_gt_node1.push_back ( sp_node.child[i]->Net_node_contains_gt_node2[j]);
                    dout << &this->NodeContainer[sp_node.Net_node_contains_gt_node1.back()] << ", ";
                }
            }
            //for ( size_t j = 0; j < sp_node.child[i]->Net_node_contains_gt_node1.size(); j++){
                //bool unique = true;
                //for (size_t check_i = 0; check_i < sp_node.Net_node_contains_gt_node1.size(); check_i++) {
                    //if (sp_node.child[i]->Net_node_contains_gt_node1[j] == sp_node.Net_node_contains_gt_node1[check_i]){
                        //unique = false;
                        //break;
                    //}
                //}

                //if (unique){
                    //sp_node.Net_node_contains_gt_node1.push_back ( sp_node.child[i]->Net_node_contains_gt_node1[j]);
                    //dout << &this->NodeContainer[sp_node.Net_node_contains_gt_node1.back()] << ", ";
                //}
            //}
        }
        else continue;
    }
    dout << sp_node.Net_node_contains_gt_node1.size() << " lineages in total" << endl;
}


void simTree::assign_lineages_at_sp_node(Node &sp_node){
    dout<<"hybrid node"<<endl;
    vector < size_t > index_container_tmp = sp_node.Net_node_contains_gt_node1;
    sp_node.Net_node_contains_gt_node1.clear();
    double left_para = sp_node.extract_hybrid_para();
    for ( size_t i = 0; i < index_container_tmp.size(); i++ ){
        if ( unifRand() < left_para ) {sp_node.Net_node_contains_gt_node1.push_back(index_container_tmp[i]);}
        else                          {sp_node.Net_node_contains_gt_node2.push_back(index_container_tmp[i]);}
    }
}


void simTree::compute_bl_extension( double multi_merge_para, size_t num_lineage ){
    current_N_lineage_To_Coalesce = 0; // reset
    current_lineage_Extension = 0; // reset
    if ( num_lineage <= 1) return;

    if ( multi_merge_para != 2 ) this->build_nc_X( num_lineage );
    current_N_lineage_To_Coalesce = (multi_merge_para == 2) ? 2 : this->update_nc();
    current_lineage_Extension = (multi_merge_para == 2) ? kingman_bl(num_lineage) : this->nc_X.min();
    assert(current_N_lineage_To_Coalesce > 1);
}


/*! \brief Write a fixed parameter into a externed newick formatted network string*/
string write_para_into_tree( string in_str /*! Externed newick formatted network string*/,
                             double para /*! Coalescent parameter or fixed population sizes */){
    if ( in_str.size() == 0 ) throw std::invalid_argument("Please define the input tree (network).");

    Net para_Net( in_str );
    for ( size_t i = 0; i < para_Net.NodeContainer.size(); i++){
        para_Net.NodeContainer[i].set_brchlen1( para );
        if ( para_Net.NodeContainer[i].hybrid() )  para_Net.NodeContainer[i].set_brchlen2( para );
    }
    para_Net.rewrite_node_content();
    return construct_adding_new_Net_str(para_Net);
}


string construct_adding_new_Net_str( Net &in_Net ){
    string out_str = in_Net.NodeContainer.back().node_content + in_Net.NodeContainer.back().label;
    if ( in_Net.NodeContainer.back().brchlen1() != 0 ) {
        out_str += ":" + to_string( static_cast<double>(in_Net.NodeContainer.back().brchlen1() ));
    }
    out_str.push_back(';');
    return out_str;
}
