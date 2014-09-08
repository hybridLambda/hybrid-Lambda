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

#include"sim_gt.hpp"
#include"mtrand.h"

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
}


void SimulationParameters::finalize(){
    // If coalescent parameter is ungiven, use Kingman coalescent as default
    if ( !this->mm_bool ) {
        para_string = write_para_into_tree( net_str, 2.0 ); 
        clog << "Default Kingman coalescent on all branches." << endl;
    }
    // If the population size is ungiven, use default population size 10000
	if ( !this->pop_bool ) {
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
    
    this->total_num_lineage = 0;
    for ( size_t i = 0; i < this->sample_size.size();i++){
		this->total_num_lineage += this->sample_size[i];
	}
    
	sp_string_pop_size = rewrite_pop_string_by_para_string(para_string,sp_string_pop_size);  // checking if modify pop_size_string is needed,

	if ( this->num_gener_bool ) net_str = write_sp_string_in_coal_unit(net_str, sp_string_pop_size);	// Convert number of generations and population size to coalescent unit
	
    sp_string_coal_unit = net_str;

    }


/*! \brief Beta function, requires tgamma function from math.h \return double */
double Beta(double x,double y){
	double Beta_return;
//	Beta_return=tgamma(x)*tgamma(y)/tgamma(x+y);
    Beta_return=exp(log(tgamma(x))+log(tgamma(y))-log(tgamma(x+y)));
	return Beta_return;
}


/*! \fn double unifRand()
 * \brief Simulate random variable between 0 and 1.
 */
double unifRand(){
	MTRand_closed return_value;
	return return_value();
    //return rand()/(double(RAND_MAX)+1);
    //return rand() / double(RAND_MAX); //generates a psuedo-random float between 0.0 and 0.999...
} 


/*! \fn int poisson_rand_var(double lambda)
 * \brief Simulating Poisson random variable from given lambda
 */
int poisson_rand_var(double lambda){
	double L=exp(-lambda);
	int k=0;
	double p=1;
	while (p>L){
         k =k + 1;
         p=p*unifRand();
	}
	//Generate uniform random number u in [0,1] and let pp × u.
    k=k-1;	
	return k;
}


sim_one_gt::sim_one_gt ( SimulationParameters* sim_param, action_board* simulation_jobs, ofstream &Si_table): parameters_(sim_param), simulation_jobs_ (simulation_jobs) {
	dout<<"	Starting simulating gene tree from "<<  this->parameters_->sp_string_coal_unit<<endl;
	//dout<<"	begin of sim_one_gt::sim_one_gt(sim::param sim_param,action_board my_action)"<<endl;
    this->Si_table_ = &Si_table;
	string sp_string_coal_unit= this->parameters_->sp_string_coal_unit;
	string sp_string_pop_size= this->parameters_->sp_string_pop_size;
	string para_string= this->parameters_->para_string;
	vector < int > sample_size= this->parameters_->sample_size;
	double mutation_rate= this->parameters_->mutation_rate;
	
	bool sim_num_gener_bool= this->simulation_jobs_->sim_num_gener_bool;
	if ( this->simulation_jobs_->sim_mut_unit_bool ||  this->simulation_jobs_->sim_num_mut_bool){
		sim_num_gener_bool=true;
	}
	Net my_Net(sp_string_coal_unit);
	Net my_pop_net(sp_string_pop_size);
	Net my_para_net(para_string);
	//vector <Node*> Net_node_ptr;
	//for ( size_t node_i=0;node_i<my_Net.NodeContainer.size();node_i++){
		//Node* new_node_ptr=NULL;
        //Net_node_ptr.push_back(new_node_ptr);
        //Net_node_ptr[node_i]=&my_Net.NodeContainer[node_i];
	//}
	tax_name=my_Net.tax_name;
	Net my_gt_coal_unit;
	my_gt_coal_unit.is_Net=false;
	
	Net my_gt_num_gener;
	if (sim_num_gener_bool){
		my_gt_num_gener.is_Net=false;
	}
	
	for ( size_t i=0;i < my_Net.NodeContainer.size();i++){
		if (my_Net.NodeContainer[i].tip_bool){
			for ( size_t sample_size_i=0;sample_size_i<sample_size.size();sample_size_i++){
				if (my_Net.descndnt[i][sample_size_i]>0){
					for (int pop_i=0;pop_i<sample_size[sample_size_i];pop_i++){
						my_gt_coal_unit.NodeContainer.push_back(my_Net.NodeContainer[i]);
						my_gt_coal_unit.descndnt.push_back(my_Net.descndnt[i]);
						ostringstream pop_i_size_str;
						pop_i_size_str<<pop_i+1;
						my_gt_coal_unit.NodeContainer.back().label=my_gt_coal_unit.NodeContainer.back().label+"_"+pop_i_size_str.str();
						my_gt_coal_unit.NodeContainer.back().brchlen1=0;
						my_gt_coal_unit.NodeContainer.back().parent1=NULL;
						my_gt_coal_unit.NodeContainer.back().parent2=NULL;
						if (sim_num_gener_bool){
							my_gt_num_gener.NodeContainer.push_back(my_gt_coal_unit.NodeContainer.back());
							my_gt_num_gener.descndnt.push_back(my_gt_coal_unit.descndnt.back());
						}
					}
				}
			}
		}	
	}
	
	vector < size_t> remaining_sp_node;
	for ( size_t sp_node_i=0;sp_node_i < my_Net.NodeContainer.size();sp_node_i++){
		if (!my_Net.NodeContainer[sp_node_i].tip_bool){
			remaining_sp_node.push_back(sp_node_i);
		}
		else{
			if (sample_size.size()>1){
				remaining_sp_node.push_back(sp_node_i);
			}
			for ( size_t gt_node_i=0;gt_node_i<my_gt_coal_unit.NodeContainer.size();gt_node_i++){
				valarray <bool> comp = (my_gt_coal_unit.descndnt[gt_node_i]==my_Net.descndnt[sp_node_i]);
				if ( comp.min() == true ){
					my_Net.NodeContainer[sp_node_i].Net_node_contains_gt_node1.push_back(gt_node_i);
				}
			}
		}
	}
	int num_tax=my_Net.tax_name.size();
	vector < size_t> remaining_gt_node;
	int gt_num_tips=my_gt_coal_unit.NodeContainer.size();
	for (int interior_i=1;interior_i<gt_num_tips;interior_i++){
		Node new_interior_node;
		string interior_node_label_str;
		if (interior_i==(gt_num_tips-1)){
			interior_node_label_str="root";
		}
		else{
			interior_node_label_str="Int_";
			ostringstream interior_i_str;
			interior_i_str<<interior_i;
			interior_node_label_str=interior_node_label_str+interior_i_str.str();
		}
		new_interior_node.label=interior_node_label_str;
		remaining_gt_node.push_back(my_gt_coal_unit.NodeContainer.size());//my_gt_coal_unit.NodeContainer.size() is the current index of the gt interior node.
		my_gt_coal_unit.NodeContainer.push_back(new_interior_node);
		valarray <int> intialize_descndnt(num_tax);
		my_gt_coal_unit.descndnt.push_back(intialize_descndnt);
		if (sim_num_gener_bool){
			my_gt_num_gener.NodeContainer.push_back(my_gt_coal_unit.NodeContainer.back());
			my_gt_num_gener.descndnt.push_back(my_gt_coal_unit.descndnt.back());
		}
	}
	
	vector <Node*> gt_nodes_ptr;
	for ( size_t gt_node_i=0;gt_node_i<my_gt_coal_unit.NodeContainer.size();gt_node_i++){
		Node* new_node_ptr=NULL;
		gt_nodes_ptr.push_back(new_node_ptr);
		gt_nodes_ptr[gt_node_i]=&my_gt_coal_unit.NodeContainer[gt_node_i];
		valarray <int> descndnt2_dummy;
		my_gt_coal_unit.descndnt2.push_back(descndnt2_dummy);
	}
	vector <Node*> num_gener_gt_nodes_ptr;
	if (sim_num_gener_bool){
		for ( size_t gt_node_i=0;gt_node_i<my_gt_num_gener.NodeContainer.size();gt_node_i++){
			Node* new_node_ptr=NULL;
			num_gener_gt_nodes_ptr.push_back(new_node_ptr);
			num_gener_gt_nodes_ptr[gt_node_i]=&my_gt_num_gener.NodeContainer[gt_node_i];
			valarray <int> descndnt2_dummy;
			my_gt_num_gener.descndnt2.push_back(descndnt2_dummy);
	
		}
	}
	int rank_i=1;

	 size_t remaining_sp_node_i=0;	
	while (remaining_sp_node.size()>0){
		
		int node_i=remaining_sp_node[remaining_sp_node_i];
		if (my_Net.NodeContainer[node_i].rank() == rank_i){
			for ( size_t child_i=0;child_i<my_Net.NodeContainer[node_i].child.size();child_i++){
				if (my_Net.NodeContainer[node_i].child[child_i]->parent1->label==my_Net.NodeContainer[node_i].label){
					for ( size_t child_i_contains_gt_node1_i=0;child_i_contains_gt_node1_i<my_Net.NodeContainer[node_i].child[child_i]->Net_node_contains_gt_node1.size();child_i_contains_gt_node1_i++){
						my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.push_back(my_Net.NodeContainer[node_i].child[child_i]->Net_node_contains_gt_node1[child_i_contains_gt_node1_i]);
					}
				}
				else{
					if (my_Net.NodeContainer[node_i].child[child_i]->parent2->label!=my_Net.NodeContainer[node_i].label){
						throw std::invalid_argument("parent2 is not correct!!!");
					}
					for ( size_t child_i_contains_gt_node1_i=0;child_i_contains_gt_node1_i<my_Net.NodeContainer[node_i].child[child_i]->Net_node_contains_gt_node2.size();child_i_contains_gt_node1_i++){
						my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.push_back(my_Net.NodeContainer[node_i].child[child_i]->Net_node_contains_gt_node2[child_i_contains_gt_node1_i]);
					}					
				}
			}

			dout<<"In population: "<< my_Net.NodeContainer[node_i].label<<" with rank  "<<rank_i;
			
			if (rank_i==my_Net.max_rank){dout<<", at the root, everything coalesces"<<endl;}else{dout<<endl;}
			//here to choose go left or right for hybrid node.			
			if (my_Net.NodeContainer[node_i].parent2){
				dout<<"hybrid node"<<endl;
				vector < size_t > Net_node_contains_gt_node_dummy=my_Net.NodeContainer[node_i].Net_node_contains_gt_node1;
				my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.clear();
				double left_para=extract_hybrid_para(my_Net.NodeContainer[node_i].label);
				for ( size_t contains_gt_node_i=0;contains_gt_node_i<Net_node_contains_gt_node_dummy.size();contains_gt_node_i++){
					if (unifRand()<left_para){
						my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.push_back(Net_node_contains_gt_node_dummy[contains_gt_node_i]);
					}
					else{
						my_Net.NodeContainer[node_i].Net_node_contains_gt_node2.push_back(Net_node_contains_gt_node_dummy[contains_gt_node_i]);
					}
				}
			}
			
			double num_lineage=double(my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.size());
						
			vector < vector <double> > lambda_bk_mat;
			if (my_para_net.NodeContainer[node_i].brchlen1!=2){
				lambda_bk_mat=build_lambda_bk_mat(my_para_net.NodeContainer[node_i].brchlen1,num_lineage);			
			}
			
			double remaining_length=my_Net.NodeContainer[node_i].brchlen1;
			dout<<num_lineage<<" lineages entering population "<< my_Net.NodeContainer[node_i].label <<" with remaining branch length of "<<remaining_length  <<endl;

			int nc=0;
			double length=0;
			if (num_lineage>1 && nc <2){
				if (my_para_net.NodeContainer[node_i].brchlen1!=2){
					valarray <double> nc_X=build_nc_X(lambda_bk_mat, num_lineage);
					nc=update_nc(nc_X);
					length=nc_X.min();											
				}
				else{
					nc=2;
					length=kingman_bl(num_lineage);
				}
			}
			
			//dout<<nc<<" lineages coalesce "<<endl;
				//dout<<"(0)num_lineage "<<num_lineage<<endl;
				//dout<<"(0)remaining length is   "<<remaining_length<<endl;
				//dout<<"(0)proposed exp length is   "<<length<<endl;
			int lineage_index;
			 size_t gt_child_node_index;
			while (((length < remaining_length) && (num_lineage>1) )  || ( (rank_i==my_Net.max_rank) && (num_lineage>1 )) ){
					
					//dout<<endl<<endl;
					//dout<<"(1)num_lineage "<<num_lineage<<endl;
					//dout<<"(1)remaining length is   "<<remaining_length<<endl;
					//dout<<"(1)proposed exp length is   "<<length<<endl;		
				if (nc>1){
					dout<<"  "<<nc<<" lineages coalesce at time "<<length<<endl;
					for (int nc_i=0;nc_i<nc;nc_i++){
						lineage_index= rand() % my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.size();
						gt_child_node_index=my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[lineage_index];
						gt_nodes_ptr[gt_child_node_index]->brchlen1=gt_nodes_ptr[gt_child_node_index]->brchlen1+length; // b_alpha <- b_alpha + l

						add_node(gt_nodes_ptr[remaining_gt_node[0]],gt_nodes_ptr[gt_child_node_index]);
						if (sim_num_gener_bool){
							num_gener_gt_nodes_ptr[gt_child_node_index]->brchlen1=num_gener_gt_nodes_ptr[gt_child_node_index]->brchlen1+(length * my_pop_net.NodeContainer[node_i].brchlen1) ; // b_alpha <- b_alpha + l*pop_size
							add_node(num_gener_gt_nodes_ptr[remaining_gt_node[0]],num_gener_gt_nodes_ptr[gt_child_node_index]);
						}
						my_gt_coal_unit.descndnt[remaining_gt_node[0]]=my_gt_coal_unit.descndnt[remaining_gt_node[0]]+my_gt_coal_unit.descndnt[gt_child_node_index];
						if (sim_num_gener_bool){
							my_gt_num_gener.descndnt[remaining_gt_node[0]]=my_gt_num_gener.descndnt[remaining_gt_node[0]]+my_gt_num_gener.descndnt[gt_child_node_index];
						}
						my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.erase(my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.begin()+lineage_index); // X'\{alpha}
					}
				}
				
				my_gt_coal_unit.NodeContainer[remaining_gt_node[0]].num_descndnt=my_gt_coal_unit.descndnt[remaining_gt_node[0]].sum();
				if (sim_num_gener_bool){
					my_gt_num_gener.NodeContainer[remaining_gt_node[0]].num_descndnt=my_gt_num_gener.descndnt[remaining_gt_node[0]].sum();
				}
				my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.push_back(remaining_gt_node[0]);// introducing a new lineage gamma


				gt_nodes_ptr[remaining_gt_node[0]]->height=gt_nodes_ptr[remaining_gt_node[0]]->child[0]->height+gt_nodes_ptr[remaining_gt_node[0]]->child[0]->brchlen1; // a_gamma <- a_alpha + b_alpha
				for ( size_t update_brch_i=0;update_brch_i<my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.size()-1;update_brch_i++){
					gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[update_brch_i]]->brchlen1=gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[update_brch_i]]->brchlen1+length;		
					if (sim_num_gener_bool){
						num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[update_brch_i]]->brchlen1= length*my_pop_net.NodeContainer[node_i].brchlen1+num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[update_brch_i]]->brchlen1;
					}
				}
				remaining_gt_node.erase(remaining_gt_node.begin());
				
				remaining_length= remaining_length - length;
					//dout<<"remaining_length "<<remaining_length<<endl;
				num_lineage=double(my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.size());

					//dout<<"checking num_lineage "<<num_lineage<<endl;

				nc=0;
				if (num_lineage>1 && nc <2){
					if (my_para_net.NodeContainer[node_i].brchlen1!=2){
						valarray <double> nc_X=build_nc_X(lambda_bk_mat, num_lineage);
						nc=update_nc(nc_X);
						length=nc_X.min();
							//dout<<"number of lineages gonna coalesce "<<nc<<endl;
					}
					else{
						nc=2;
						length=kingman_bl(num_lineage);
					}
				}
				dout<<num_lineage<<" live lineages in population "<< my_Net.NodeContainer[node_i].label <<" with remaining branch length of "<<remaining_length  <<endl;
			
			}
			
				
			
			if (my_Net.NodeContainer[node_i].parent2){			
				dout<<"hybrid node parent1 finished, in parent 2 now"<<endl;
				double num_lineage=(my_Net.NodeContainer[node_i].Net_node_contains_gt_node2.size());
				lambda_bk_mat.clear();
				if (my_para_net.NodeContainer[node_i].brchlen2!=2){
					lambda_bk_mat=build_lambda_bk_mat(my_para_net.NodeContainer[node_i].brchlen2,num_lineage);
				}
											
				double remaining_length=my_Net.NodeContainer[node_i].brchlen2;
				dout<<num_lineage<<" lineages entering population "<< my_Net.NodeContainer[node_i].label <<" with remaining branch length of "<<remaining_length  <<endl;

				double length=0;		
				int nc=0;
				if (num_lineage>1 && nc <2){
					if (my_para_net.NodeContainer[node_i].brchlen2!=2){
						valarray <double> nc_X=build_nc_X(lambda_bk_mat, num_lineage);					
						nc=update_nc(nc_X);
						length=nc_X.min();
							//dout<<"number of lineages gonna coalesce "<<nc<<endl;		
					}
					else{
						nc=2;
						length=kingman_bl(num_lineage);
					}
				}
				
					//dout<<"length  "<<length<<endl;
					//dout<<"(0)num_lineage "<<num_lineage<<endl;
					//dout<<"(0)remaining length is   "<<remaining_length<<endl;
					//dout<<"(0)proposed exp length is   "<<length<<endl;

				int lineage_index;
				 size_t gt_child_node_index;
				while ( (length < remaining_length && num_lineage>1.0 )  ){
						//dout<<endl<<endl;
						//dout<<"(1)num_lineage "<<num_lineage<<endl;
						//dout<<"(1)remaining length is   "<<remaining_length<<endl;
						//dout<<"(1)proposed exp length is   "<<length<<endl;		

					if (nc>1){
						dout<<"  "<<nc<<" lineages coalesce at time "<<length<<endl;

						for (int nc_i=0;nc_i<nc;nc_i++){
							lineage_index= rand() % my_Net.NodeContainer[node_i].Net_node_contains_gt_node2.size();
							gt_child_node_index=my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[lineage_index];
							gt_nodes_ptr[gt_child_node_index]->brchlen1=gt_nodes_ptr[gt_child_node_index]->brchlen1+length; // b_alpha <- b_alpha + l
							add_node(gt_nodes_ptr[remaining_gt_node[0]],gt_nodes_ptr[gt_child_node_index]);
							if (sim_num_gener_bool){
								num_gener_gt_nodes_ptr[gt_child_node_index]->brchlen1=num_gener_gt_nodes_ptr[gt_child_node_index]->brchlen1+length*my_pop_net.NodeContainer[node_i].brchlen2; // b_alpha <- b_alpha + l
								add_node(num_gener_gt_nodes_ptr[remaining_gt_node[0]],num_gener_gt_nodes_ptr[gt_child_node_index]);
							}
							my_gt_coal_unit.descndnt[remaining_gt_node[0]]=my_gt_coal_unit.descndnt[remaining_gt_node[0]]+my_gt_coal_unit.descndnt[gt_child_node_index];
							if (sim_num_gener_bool){
								my_gt_num_gener.descndnt[remaining_gt_node[0]]=my_gt_num_gener.descndnt[remaining_gt_node[0]]+my_gt_num_gener.descndnt[gt_child_node_index];
							}
							my_Net.NodeContainer[node_i].Net_node_contains_gt_node2.erase(my_Net.NodeContainer[node_i].Net_node_contains_gt_node2.begin()+lineage_index); // X'\{alpha}
							
						}
					}

					my_gt_coal_unit.NodeContainer[remaining_gt_node[0]].num_descndnt=my_gt_coal_unit.descndnt[remaining_gt_node[0]].sum();
					if (sim_num_gener_bool){
						my_gt_num_gener.NodeContainer[remaining_gt_node[0]].num_descndnt=my_gt_num_gener.descndnt[remaining_gt_node[0]].sum();
					}
					my_Net.NodeContainer[node_i].Net_node_contains_gt_node2.push_back(remaining_gt_node[0]);

					gt_nodes_ptr[remaining_gt_node[0]]->height=gt_nodes_ptr[remaining_gt_node[0]]->child[0]->height+gt_nodes_ptr[remaining_gt_node[0]]->child[0]->brchlen1;
					for ( size_t update_brch_i=0;update_brch_i<my_Net.NodeContainer[node_i].Net_node_contains_gt_node2.size()-1;update_brch_i++){
						gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[update_brch_i]]->brchlen1=gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[update_brch_i]]->brchlen1+length;		
						if (sim_num_gener_bool){
							num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[update_brch_i]]->brchlen1= length*my_pop_net.NodeContainer[node_i].brchlen2 + num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[update_brch_i]]->brchlen1;
						}
					}
					remaining_gt_node.erase(remaining_gt_node.begin());
					remaining_length= remaining_length - length;
					num_lineage=double(my_Net.NodeContainer[node_i].Net_node_contains_gt_node2.size());
						//dout<<"checking num_lineage "<<num_lineage<<endl;
					
					nc=0;
					if (num_lineage>1 && nc <2){
						if (my_para_net.NodeContainer[node_i].brchlen2!=2){
							valarray <double> nc_X=build_nc_X(lambda_bk_mat, num_lineage);						
							nc=update_nc(nc_X);
							length=nc_X.min();
						}
						else{
							nc=2;
							length=kingman_bl(num_lineage);
						}
					}

						//dout<<"length  "<<length<<endl;
						//dout<<"(0)num_lineage "<<num_lineage<<endl;
						//dout<<"(0)remaining length is   "<<remaining_length<<endl;
						//dout<<"(0)proposed exp length is   "<<length<<endl;
					dout<<num_lineage<<" live lineages in population "<< my_Net.NodeContainer[node_i].label <<" with remaining branch length of "<<remaining_length  <<endl;

				}
			}
					
			dout<<"*************************before adjusting***************"<<endl;
			assert(my_gt_coal_unit.print_all_node_dout());
			
			if (rank_i<my_Net.max_rank){
				if (my_Net.NodeContainer[node_i].parent2){
					for ( size_t num_lineage_i=0;num_lineage_i<my_Net.NodeContainer[node_i].Net_node_contains_gt_node2.size();num_lineage_i++){
					//	gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node[num_lineage_i]]->brchlen1=gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node[num_lineage_i]]->brchlen1+my_Net.NodeContainer[node_i].parent1->height-gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node[num_lineage_i]]->height;
					//	gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node[num_lineage_i]]->height=my_Net.NodeContainer[node_i].parent1->height;	
						if (sim_num_gener_bool){
							if ((gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[num_lineage_i]]->height)>  my_Net.NodeContainer[node_i].height){
							  num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[num_lineage_i]]->brchlen1= (my_Net.NodeContainer[node_i].parent2->height - gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[num_lineage_i]]->height)*my_pop_net.NodeContainer[node_i].brchlen2;//+num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1;						
							}
							else{
							  num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[num_lineage_i]]->brchlen1= (my_Net.NodeContainer[node_i].parent2->height - gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[num_lineage_i]]->height- gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[num_lineage_i]]->brchlen1)*my_pop_net.NodeContainer[node_i].brchlen2+num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[num_lineage_i]]->brchlen1;						
							}
						}
						gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[num_lineage_i]]->brchlen1=my_Net.NodeContainer[node_i].parent2->height -gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node2[num_lineage_i]]->height;
					}					
				}
				//else{
				// through parent 1, it should always be checked!!!	
				
				for ( size_t num_lineage_i=0;num_lineage_i<my_Net.NodeContainer[node_i].Net_node_contains_gt_node1.size();num_lineage_i++){
					if (sim_num_gener_bool){
						if ((gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[num_lineage_i]]->height)>  my_Net.NodeContainer[node_i].height){
							num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1= (my_Net.NodeContainer[node_i].parent1->height - gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[num_lineage_i]]->height)*my_pop_net.NodeContainer[node_i].brchlen1;//+num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1;						
						}
						else{
							num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1= (my_Net.NodeContainer[node_i].parent1->height - gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[num_lineage_i]]->height- gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1)*my_pop_net.NodeContainer[node_i].brchlen1+num_gener_gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1;						
						}
					}
					gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1=my_Net.NodeContainer[node_i].parent1->height -gt_nodes_ptr[my_Net.NodeContainer[node_i].Net_node_contains_gt_node1[num_lineage_i]]->height;
				}
				//}
				
					dout<<"************************* after adjusting***************"<<endl;
					assert(my_gt_coal_unit.print_all_node_dout());

			}
			remaining_sp_node.erase(remaining_sp_node.begin()+remaining_sp_node_i);
				//dout<<rank_i<<" "<<my_Net.NodeContainer[node_i].label<<" "<<my_Net.NodeContainer[node_i].path_time.size()<<endl;
		}
		else{
			remaining_sp_node_i++;
		}

		if (remaining_sp_node_i==remaining_sp_node.size()-1){
			rank_i++;
			remaining_sp_node_i=0;
		}
	}
	
	for ( size_t node_i=0;node_i<gt_nodes_ptr.size();){
		if (gt_nodes_ptr[node_i]->num_descndnt==0){
			gt_nodes_ptr.erase(gt_nodes_ptr.begin()+node_i);
			if (sim_num_gener_bool){
				num_gener_gt_nodes_ptr.erase(num_gener_gt_nodes_ptr.begin()+node_i);
			}
		}
		else{
			node_i++;
		}
	}
	
	
	for ( size_t node_i=0;node_i<my_gt_coal_unit.NodeContainer.size();node_i++){
		if (my_gt_coal_unit.descndnt[node_i].sum()==0){
			break;
		}
		if (!my_gt_coal_unit.NodeContainer[node_i].tip_bool){
			for ( size_t tax_name_i=0;tax_name_i<tax_name.size();tax_name_i++){
				if (my_gt_coal_unit.descndnt[node_i][tax_name_i] > 0 ){
					if (my_gt_coal_unit.NodeContainer[node_i].clade.size()==0){
						my_gt_coal_unit.NodeContainer[node_i].clade=tax_name[tax_name_i];
					}
					else{
						my_gt_coal_unit.NodeContainer[node_i].clade=my_gt_coal_unit.NodeContainer[node_i].clade+tax_name[tax_name_i];
					}
					my_gt_coal_unit.NodeContainer[node_i].clade.push_back('&');
				}
			}
			my_gt_coal_unit.NodeContainer[node_i].clade.erase(my_gt_coal_unit.NodeContainer[node_i].clade.size()-1,1);
		}
	}
	

	//if (!my_gt_coal_unit.is_ultrametric_func()){throw "not ultrametic";}
	

	rewrite_node_content(gt_nodes_ptr);
	string old_gt_string_coal_unit=gt_nodes_ptr.back()->node_content+gt_nodes_ptr.back()->label+";";

	if (sim_num_gener_bool){
		rewrite_node_content(num_gener_gt_nodes_ptr);
		string old_gt_string_num_gener=num_gener_gt_nodes_ptr.back()->node_content+num_gener_gt_nodes_ptr.back()->label+";";
		gt_string_gener_num=remove_interior_label(old_gt_string_num_gener);
	}

	//dout<<"flag 2"<<endl;
	//dout<<gt_nodes_ptr.back()->label<<endl;
	//if (sim_num_gener_bool){
		//dout<<num_gener_gt_nodes_ptr.back()->label<<endl;
	//}
	//dout<<old_gt_string_coal_unit<<endl;

	gt_string_coal_unit=remove_interior_label(old_gt_string_coal_unit);
	// check if the gene tree is ultramatric.
	dout<<"check of if "<<gt_string_coal_unit <<" is ultrametric"<<endl;
	Net checking_ultra_net(gt_string_coal_unit);
	//if (!checking_ultra_net.is_ultrametric){throw "Gene tree is not ultrametric";}
	
	
	if ( this->simulation_jobs_->sim_mut_unit_bool){
		build_gt_string_mut_unit(mutation_rate);
	}	
	
	
	if ( this->simulation_jobs_->sim_num_mut_bool ||  this->simulation_jobs_->Si_num_bool){
		Net mt_tree(gt_string_gener_num);
		vector <double> brch_total;
		vector <Node*> mt_nodes_ptr;
		total_brchlen=0;
		for ( size_t gt_node_i=0;gt_node_i<mt_tree.NodeContainer.size();gt_node_i++){
			Node* new_node_ptr=NULL;
			mt_nodes_ptr.push_back(new_node_ptr);
			total_brchlen=total_brchlen+mt_tree.NodeContainer[gt_node_i].brchlen1;
			mt_nodes_ptr[gt_node_i]=&mt_tree.NodeContainer[gt_node_i];
			mt_nodes_ptr[gt_node_i]->brchlen1=0;
			brch_total.push_back(total_brchlen);
		}	
		
		brch_total.back()=0;
		//cout<<total_brchlen<<endl;
		//double theta=1.0;
		this->total_mut = poisson_rand_var(mutation_rate*total_brchlen);
		for ( int mut_i=0; mut_i < this->total_mut; mut_i++){
			 size_t brch_index=0;
			double u=unifRand()*total_brchlen;
			while (u>brch_total[brch_index]){
				 brch_index =brch_index + 1;	
			}
			mt_nodes_ptr[brch_index]->brchlen1 = mt_nodes_ptr[brch_index]->brchlen1+1;
		}
		rewrite_node_content(mt_nodes_ptr);
		gt_string_mut_num=mt_nodes_ptr.back()->node_content+mt_nodes_ptr.back()->label+";";
		gt_string_mut_num=remove_interior_label(gt_string_mut_num);
		
		if ( this->simulation_jobs_->Si_num_bool){	
			Si_num_out_table(mt_tree);
		}
	}

	if (tax_name.size()==2 &&  this->simulation_jobs_->mono_bool){
		compute_monophyly_vec(my_gt_coal_unit,sample_size);
	}
	dout<<"	end of sim_one_gt::sim_one_gt(sim::param sim_param,action_board my_action)"<<endl;
}


void sim_one_gt::Si_num_out_table(Net mt_tree){
	vector <int> S_i(mt_tree.tip_name.size()-1,0);
	for ( size_t sii=0;sii<S_i.size();sii++){
		for ( size_t node_i=0;node_i<mt_tree.NodeContainer.size();node_i++){
			int sii_pluse_1=sii+1;
			if ((mt_tree.NodeContainer[node_i].brchlen1>0) && (sii_pluse_1==mt_tree.descndnt2[node_i].sum() )){
				//if ((sii+1)==mt_tree.descndnt2[node_i].sum()){
					S_i[sii]=S_i[sii]+mt_tree.NodeContainer[node_i].brchlen1;
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


vector < vector <double> > build_lambda_bk_mat(double para,double num_lineage){
	vector < vector <double> > lambda_bk_mat;
	for (double b_i=2;b_i<=num_lineage;b_i++){
		vector <double> lambda_bk_mat_b;
		for (double k_i=2;k_i<=b_i;k_i++){
			double lambda_bk_mat_b_k;
			if (para<1){//0<psi<1
				//lambda_bk_mat_b_k=n_choose_k(b_i,k_i)*pow(para,k_i)*pow(1-para,b_i-k_i);//.2 is psi lambda_bk=\binom{b}{k}\psi^k (1-\psi)^{b-k}
                //lambda_bk_mat_b_k=exp(log(n_choose_k(b_i,k_i)) + log(pow(para,k_i)) + log(pow(1-para,b_i-k_i)));//.2 is psi lambda_bk=\binom{b}{k}\psi^k (1-\psi)^{b-k}
                //cout<<"normal calculation : "<<exp(log(n_choose_k(b_i,k_i)) + log(pow(para,k_i)) + log(pow(1-para,b_i-k_i)))<<endl;
                lambda_bk_mat_b_k=exp(log(boost::math::binomial_coefficient<double>(unsigned(b_i),unsigned(k_i))) + log(pow(para,k_i)) + log(pow(1-para,b_i-k_i)));//.2 is psi lambda_bk=\binom{b}{k}\psi^k (1-\psi)^{b-k}
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
                lambda_bk_mat_b_k=exp(log(boost::math::binomial_coefficient<double>(unsigned(b_i),unsigned(k_i)))+log(Beta(k_i-para,b_i-k_i+para)) - log(Beta(2.0-para,para)));// \lambda_{bk}=\binom{b}{k}\frac{B(k-\alpha,b-k+\alpha)}{B(2-\alpha,\alpha)}
				//cout<<"boost result "<<lambda_bk_mat_b_k<<endl;
			}
			lambda_bk_mat_b.push_back(lambda_bk_mat_b_k);
			}
		lambda_bk_mat.push_back(lambda_bk_mat_b);
	}
	return lambda_bk_mat;
}

void sim_one_gt::build_gt_string_mut_unit(double mutation_rate){
	Net gt_mut_unit(gt_string_gener_num);
	vector <Node*> gt_mut_unit_nodes_ptr;
	for ( size_t gt_node_i=0;gt_node_i<gt_mut_unit.NodeContainer.size();gt_node_i++){
		Node* new_node_ptr=NULL;
		gt_mut_unit_nodes_ptr.push_back(new_node_ptr);
		gt_mut_unit_nodes_ptr[gt_node_i]=&gt_mut_unit.NodeContainer[gt_node_i];
		gt_mut_unit_nodes_ptr[gt_node_i]->brchlen1=gt_mut_unit_nodes_ptr[gt_node_i]->brchlen1 * mutation_rate;
	}
	rewrite_node_content(gt_mut_unit_nodes_ptr);
	gt_string_mut_unit=gt_mut_unit_nodes_ptr.back()->node_content+gt_mut_unit_nodes_ptr.back()->label+";";
	gt_string_mut_unit=remove_interior_label(gt_string_mut_unit);
}


void sim_one_gt::compute_monophyly_vec(Net my_gt_coal_unit,vector < int > sample_size){
	vector <double> monophyly_initial(6,0);
	monophyly = monophyly_initial;
	for ( size_t tax_i=0;tax_i<2;tax_i++){
		for ( size_t node_i=0;node_i<my_gt_coal_unit.NodeContainer.size();node_i++){
			if (my_gt_coal_unit.descndnt[node_i].sum()==0){
				break;
			}
			if ((my_gt_coal_unit.descndnt[node_i][tax_i]==sample_size[tax_i]) &&  ( my_gt_coal_unit.descndnt[node_i][tax_i]==my_gt_coal_unit.descndnt[node_i].sum() ) ){
				monophyly[tax_i]=1;
				break;
			}
		}
	}
	if ( monophyly[0] == 1 ){
		if ( monophyly[1] == 1 ){
			monophyly[2] = 1;
		}
		else{
			monophyly[3] = 1;
		}	
	}
	else{
		if (monophyly[1] == 1){
			monophyly[4] = 1;
		}
		else{
			monophyly[5] = 1;
		}
	}

	for (size_t mono_i = 0; mono_i < 6; mono_i ++){
		dout<<monophyly[mono_i]<<"  ";
	}
	dout<<endl;
}


/*! \brief Convert the network branch lengths from number of generations to coalescent unit, by dividing the population size 
 * \return string Network in extended newick form, branch lengths are in coalescent unit*/
string write_sp_string_in_coal_unit(string sp_num_gener_string /*! Network in extended newick form, branch lengths are the number of generations*/,
string pop_size_string /*! Network in extended newick form, branch lengths are the population sizes*/){
	Net sp_num_gener_net(sp_num_gener_string);
	Net pop_size_net(pop_size_string);
	
	vector <Node*> sp_num_gener_net_node_ptr;
	for ( size_t node_i=0;node_i<sp_num_gener_net.NodeContainer.size();node_i++){
		Node* new_node_ptr=NULL;
        sp_num_gener_net_node_ptr.push_back(new_node_ptr);
        sp_num_gener_net_node_ptr[node_i]=&sp_num_gener_net.NodeContainer[node_i];
		if (node_i<sp_num_gener_net.NodeContainer.size()-1){
			sp_num_gener_net_node_ptr[node_i]->brchlen1=(sp_num_gener_net_node_ptr[node_i]->brchlen1) / pop_size_net.NodeContainer[node_i].brchlen1;
			if (pop_size_net.NodeContainer[node_i].hybrid){
				sp_num_gener_net_node_ptr[node_i]->brchlen2=(sp_num_gener_net_node_ptr[node_i]->brchlen2) / pop_size_net.NodeContainer[node_i].brchlen2;
			}
		}
	}
	rewrite_node_content(sp_num_gener_net_node_ptr);
	string sp_coal_unit_string=construct_adding_new_Net_str(sp_num_gener_net);
	
	return sp_coal_unit_string;
}


/*! \brief For alpha coalescent, modify the population size according to the multi merger parameter, N=N^(alpha-1), as mutation must scale with the coalescent timescale, in this case is N^(alpha-1)
 * \return string  Network in extended newick form, branch lengths are the population sizes */
string rewrite_pop_string_by_para_string(
	string para_string /*! Network in extended newick form, branch lengths are the coalescent parameters */,
	string pop_size_string /*! Network in extended newick form, branch lengths are the population sizes*/
	){
	Net para_net_check(para_string);
	Net pop_size_check(pop_size_string);
	vector <Node*> pop_size_node_ptr;
	for ( size_t node_i=0;node_i<para_net_check.NodeContainer.size();node_i++){
		Node* new_node_ptr=NULL;
		pop_size_node_ptr.push_back(new_node_ptr);
		pop_size_node_ptr[node_i]=&pop_size_check.NodeContainer[node_i];
		if (para_net_check.NodeContainer[node_i].brchlen1<2 && para_net_check.NodeContainer[node_i].brchlen1>1){ // rescale the number of generations for alpha
			pop_size_node_ptr[node_i]->brchlen1 = pow(pop_size_node_ptr[node_i]->brchlen1,para_net_check.NodeContainer[node_i].brchlen1-1);
		}
		if (para_net_check.NodeContainer[node_i].hybrid){
			if (para_net_check.NodeContainer[node_i].brchlen2<2 && para_net_check.NodeContainer[node_i].brchlen2>1){
				pop_size_node_ptr[node_i]->brchlen2=pow(pop_size_node_ptr[node_i]->brchlen2,para_net_check.NodeContainer[node_i].brchlen2-1);
			}
		}
	}
	rewrite_node_content(pop_size_node_ptr);
	string pop_size_string_return = construct_adding_new_Net_str(pop_size_check);
	return pop_size_string_return;
}


double update_coal_para(vector < vector <double> > lambda_bk_mat, double num_lineage){
	double coal_para=0;
	for (int lambda_bk_i=0;lambda_bk_i<=num_lineage-2;lambda_bk_i++){
		coal_para=coal_para+lambda_bk_mat[num_lineage-2][lambda_bk_i];
		dout<<"lambda_bk_mat[num_lineage-2].size() "<<lambda_bk_mat[num_lineage-2].size()<<endl;
		dout<< " !!! "<<num_lineage<<"  "<< lambda_bk_i+2<<endl;
		dout<<"coal_para "<<coal_para<<endl;
	}
	return coal_para;
}

valarray <double> build_nc_X(vector < vector <double> > lambda_bk_mat, double num_lineage){
	valarray <double> nc_X(num_lineage-1);
	for ( size_t kmerge=0;kmerge<nc_X.size();kmerge++){
		nc_X[kmerge]= -log( 1-unifRand() )/ lambda_bk_mat[num_lineage-2][kmerge];
		//cout<<nc_X[kmerge]<<endl;
		//if (isnan(nc_X[kmerge])){cout<<lambda_bk_mat[num_lineage-2][kmerge]<<endl;}
	}
	return nc_X;
}

//use heap structure for this!
int update_nc(valarray <double> nc_X){	
	for (int kmerge = 0; kmerge < int(nc_X.size()); kmerge++){
		if (nc_X[kmerge] == nc_X.min()){
			//nc=kmerge+2;
			return kmerge+2;
			//break;
		}
	}
    dout << "k merger was never found ... " << endl;
}

double kingman_bl(double num_lineage){
	return -log( 1-unifRand() ) / num_lineage / (num_lineage-1) * 2.0;
}

	
/*! \brief Write a fixed parameter into a externed newick formatted network string*/
string write_para_into_tree(string in_str /*! Externed newick formatted network string*/, 
double para /*! Coalescent parameter or fixed population sizes */){
	if ( in_str.size() == 0 ){
		throw std::invalid_argument("Please define the input tree (network).");
	}
	Net para_Net(in_str);
	vector <Node*> para_Net_node_ptr;
	for (size_t node_i=0;node_i<para_Net.NodeContainer.size();node_i++){
		Node* new_node_ptr=NULL;
        para_Net_node_ptr.push_back(new_node_ptr);
        para_Net_node_ptr[node_i]=&para_Net.NodeContainer[node_i];
		para_Net_node_ptr[node_i]->brchlen1=para;
		if (para_Net.NodeContainer[node_i].hybrid){
			para_Net_node_ptr[node_i]->brchlen2=para;
		}
	}
	para_Net_node_ptr.back()->brchlen1=para;
	rewrite_node_content(para_Net_node_ptr);
	string para_string=construct_adding_new_Net_str(para_Net);
	return para_string;	
}
