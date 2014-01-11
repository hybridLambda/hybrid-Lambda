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

/*! \file sim_gt.cpp
 * \brief Core functions used to simulate gene trees for given species network or tree */


#include"sim_gt.hpp"
#include"mtrand.h"

action_board::action_board(){
	sim_mut_unit_bool=false;
	sim_num_gener_bool=false;
	sim_num_mut_bool=false;
	mono_bool=false;
	Si_num_bool=false;
	//total_brchlen_bool=false;
	gene_tree_file="GENE_TREE";
}

action_board::action_board(int argc, char *argv[]){
	sim_mut_unit_bool=false;
	sim_num_gener_bool=false;
	sim_num_mut_bool=false;
	mono_bool=false;
	Si_num_bool=false;
	//total_brchlen_bool=false;
	gene_tree_file="GENE_TREE";		
	for(int argc_i=1; argc_i < argc;argc_i++ ){
		std::string argv_i(argv[argc_i]);
		if (argv_i=="-sim_mut_unit"){
			sim_mut_unit_bool=true;
		}
		if (argv_i=="-sim_num_gener"){
			sim_num_gener_bool=true;
		}
		if (argv_i=="-sim_num_mut" || argv_i=="-seg"){
			sim_num_mut_bool=true;
		}
		if (argv_i=="-sim_Si_num"){
			Si_num_bool=true;
			sim_num_mut_bool=true;
			check_and_remove("out_table");
		}
		//if (argv_i=="-seg"){
			////sites_data_bool=true;
			//sim_num_mut_bool=true;
		//}
		if (argv_i=="-mono"){
			mono_bool=true;
		}
		if (argv_i=="-GENE" || argv_i=="-gF"){
			gene_tree_file=argv[argc_i+1];
			argc_i++;
		}
		

	}
}
sim::param::param(){
	mutation_rate=0.00005;
	num_sim_gt=1;
	//pop_size=10000;
	//pop_size_string_bool=false;
	//mm=2.0;
	pop_bool=false;
	mm_bool=false;
	
	num_gener_bool=false;
	sp_coal_unit_bool=false;
}

sim::param::param(int argc, char *argv[]){
	bool samples_bool=false;
	mutation_rate=0.00005;
	num_sim_gt=1;
	//pop_size=10000;
	//pop_size_string_bool=false;
	//mm=2.0;
	pop_bool=false;
	mm_bool=false;
	
	num_gener_bool=false;
	sp_coal_unit_bool=false;		
	int argc_i=1;

	while( argc_i < argc ){
		std::string argv_i(argv[argc_i]);
		
		if (argv_i=="-sp_coal_unit" || argv_i=="-sp_num_gener" || argv_i=="-spcu" || argv_i=="-spng"){
			if (argv_i=="-sp_num_gener" || argv_i=="-spng"){
				num_gener_bool=true;
			}
			if (argv_i=="-sp_coal_unit" || argv_i=="-spcu"){
				sp_coal_unit_bool=true;
			}
			if (sp_coal_unit_bool && num_gener_bool){
				throw std::invalid_argument("Species tree branch length should only be in Coalescent unit or number of generations. Choose either -sp_num_gener or -sp_coal_unit");
			}
			net_str=read_input_line(argv[argc_i+1]);
			argc_i++;

		}
		
		if (argv_i=="-mu"){
			read_input_to_param<double>(argv[argc_i+1],mutation_rate);
			argc_i++;
		}
				
		if (argv_i=="-num"){
			read_input_to_param<int>(argv[argc_i+1],num_sim_gt);
			argc_i++;
		}

		if (argv_i=="-mm"){
			mm_bool=true;
			para_string=read_input_para(argv[argc_i+1],net_str);
			argc_i++;
		}
			
		if (argv_i=="-pop"){
			pop_bool=true;
			sp_string_pop_size=read_input_para(argv[argc_i+1],net_str);
			argc_i++;
		}

		if (argv_i=="-S" ){
			samples_bool=true;
			for (int argc_j=argc_i+1;argc_j<argc;argc_j++){
				string s(argv[argc_j]);
				if (!isdigit(s[0])){
					break;
				}
				istringstream pop_num_str(s);
				int pop_num;
				pop_num_str>>pop_num;
				sample_size.push_back(pop_num);
			}
		}
		argc_i++;
	}

	if (!mm_bool){
		para_string=write_para_into_tree(net_str, 2.0); // If coalescent parameter is ungiven, use Kingman coalescent as default
	}

	if (!pop_bool){
		sp_string_pop_size=write_para_into_tree(net_str, 10000.0); // If the population size is ungiven, use default population size 10000
	}

	Net net_dummy(net_str);
	
	if (!samples_bool){
		for (int i=0;i<net_dummy.tax_name.size();i++){
			sample_size.push_back(1);
		}
	}
	else{//  check the number of lineages and the number of species 
		if (sample_size.size()!=net_dummy.tax_name.size()){
			throw std::invalid_argument("Numbers of samples and numbers of species not equal!!!");
		}
	}
	
	sp_string_pop_size=rewrite_pop_string_by_para_string(para_string,sp_string_pop_size);  // checking if modify pop_size_string is needed,

	if (num_gener_bool){
		net_str=write_sp_string_in_coal_unit(net_str,sp_string_pop_size);	// Convert number of generations and population size to coalescent unit
	}
	sp_string_coal_unit=net_str;
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


//sim_one_gt::sim_one_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string,vector < int > sample_size, double mutation_rate,bool sim_mut_unit_bool, bool sim_num_gener_bool,bool sim_num_mut_bool,bool mono_bool){
//sim_one_gt::sim_one_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string,vector < int > sample_size, double mutation_rate,action_board my_action){
sim_one_gt::sim_one_gt(sim::param sim_param,action_board my_action){
	dout<<"	Starting simulating gene tree from "<< sim_param.sp_string_coal_unit<<endl;
	//dout<<"	begin of sim_one_gt::sim_one_gt(sim::param sim_param,action_board my_action)"<<endl;

	string sp_string_coal_unit=sim_param.sp_string_coal_unit;
	string sp_string_pop_size=sim_param.sp_string_pop_size;
	string para_string=sim_param.para_string;
	vector < int > sample_size=sim_param.sample_size;
	double mutation_rate=sim_param.mutation_rate;
	
	bool sim_num_gener_bool=my_action.sim_num_gener_bool;
	if (my_action.sim_mut_unit_bool || my_action.sim_num_mut_bool){
		sim_num_gener_bool=true;
	}
	Net my_Net(sp_string_coal_unit);
	Net my_pop_net(sp_string_pop_size);
	Net my_para_net(para_string);
	vector <Node*> Net_node_ptr;
	for (unsigned int node_i=0;node_i<my_Net.Net_nodes.size();node_i++){
		Node* new_node_ptr=NULL;
        Net_node_ptr.push_back(new_node_ptr);
        Net_node_ptr[node_i]=&my_Net.Net_nodes[node_i];
	}
	tax_name=my_Net.tax_name;
	Net my_gt_coal_unit;
	my_gt_coal_unit.is_net=false;
	
	Net my_gt_num_gener;
	if (sim_num_gener_bool){
		my_gt_num_gener.is_net=false;
	}
	
	for (unsigned int i=0;i < my_Net.Net_nodes.size();i++){
		if (my_Net.Net_nodes[i].tip_bool){
			for (unsigned int sample_size_i=0;sample_size_i<sample_size.size();sample_size_i++){
				if (my_Net.descndnt[i][sample_size_i]>0){
					for (int pop_i=0;pop_i<sample_size[sample_size_i];pop_i++){
						my_gt_coal_unit.Net_nodes.push_back(my_Net.Net_nodes[i]);
						my_gt_coal_unit.descndnt.push_back(my_Net.descndnt[i]);
						ostringstream pop_i_size_str;
						pop_i_size_str<<pop_i+1;
						my_gt_coal_unit.Net_nodes.back().label=my_gt_coal_unit.Net_nodes.back().label+"_"+pop_i_size_str.str();
						my_gt_coal_unit.Net_nodes.back().brchlen1=0;
						my_gt_coal_unit.Net_nodes.back().parent1=NULL;
						my_gt_coal_unit.Net_nodes.back().parent2=NULL;
						if (sim_num_gener_bool){
							my_gt_num_gener.Net_nodes.push_back(my_gt_coal_unit.Net_nodes.back());
							my_gt_num_gener.descndnt.push_back(my_gt_coal_unit.descndnt.back());
						}
					}
				}
			}
		}	
	}
	
	vector <unsigned int> remaining_sp_node;
	for (unsigned int sp_node_i=0;sp_node_i < my_Net.Net_nodes.size();sp_node_i++){
		if (!my_Net.Net_nodes[sp_node_i].tip_bool){
			remaining_sp_node.push_back(sp_node_i);
		}
		else{
			if (sample_size.size()>1){
				remaining_sp_node.push_back(sp_node_i);
			}
			for (unsigned int gt_node_i=0;gt_node_i<my_gt_coal_unit.Net_nodes.size();gt_node_i++){
				valarray <bool> comp = (my_gt_coal_unit.descndnt[gt_node_i]==my_Net.descndnt[sp_node_i]);
				if ( comp.min() == true ){
					my_Net.Net_nodes[sp_node_i].Net_node_contains_gt_node1.push_back(gt_node_i);
				}
			}
		}
	}
	int num_tax=my_Net.tax_name.size();
	vector <unsigned int> remaining_gt_node;
	int gt_num_tips=my_gt_coal_unit.Net_nodes.size();
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
		remaining_gt_node.push_back(my_gt_coal_unit.Net_nodes.size());//my_gt_coal_unit.Net_nodes.size() is the current index of the gt interior node.
		my_gt_coal_unit.Net_nodes.push_back(new_interior_node);
		valarray <int> intialize_descndnt(num_tax);
		my_gt_coal_unit.descndnt.push_back(intialize_descndnt);
		if (sim_num_gener_bool){
			my_gt_num_gener.Net_nodes.push_back(my_gt_coal_unit.Net_nodes.back());
			my_gt_num_gener.descndnt.push_back(my_gt_coal_unit.descndnt.back());
		}
	}
	
	vector <Node*> gt_nodes_ptr;
	for (unsigned int gt_node_i=0;gt_node_i<my_gt_coal_unit.Net_nodes.size();gt_node_i++){
		Node* new_node_ptr=NULL;
		gt_nodes_ptr.push_back(new_node_ptr);
		gt_nodes_ptr[gt_node_i]=&my_gt_coal_unit.Net_nodes[gt_node_i];
		valarray <int> descndnt2_dummy;
		my_gt_coal_unit.descndnt2.push_back(descndnt2_dummy);
	}
	vector <Node*> num_gener_gt_nodes_ptr;
	if (sim_num_gener_bool){
		for (unsigned int gt_node_i=0;gt_node_i<my_gt_num_gener.Net_nodes.size();gt_node_i++){
			Node* new_node_ptr=NULL;
			num_gener_gt_nodes_ptr.push_back(new_node_ptr);
			num_gener_gt_nodes_ptr[gt_node_i]=&my_gt_num_gener.Net_nodes[gt_node_i];
			valarray <int> descndnt2_dummy;
			my_gt_num_gener.descndnt2.push_back(descndnt2_dummy);
	
		}
	}
	int rank_i=1;

	unsigned int remaining_sp_node_i=0;	
	while (remaining_sp_node.size()>0){
		
		int node_i=remaining_sp_node[remaining_sp_node_i];
		if (my_Net.Net_nodes[node_i].rank==rank_i){
			for (unsigned int child_i=0;child_i<my_Net.Net_nodes[node_i].child.size();child_i++){
				if (my_Net.Net_nodes[node_i].child[child_i]->parent1->label==my_Net.Net_nodes[node_i].label){
					for (unsigned int child_i_contains_gt_node1_i=0;child_i_contains_gt_node1_i<my_Net.Net_nodes[node_i].child[child_i]->Net_node_contains_gt_node1.size();child_i_contains_gt_node1_i++){
						my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.push_back(my_Net.Net_nodes[node_i].child[child_i]->Net_node_contains_gt_node1[child_i_contains_gt_node1_i]);
					}
				}
				else{
					if (my_Net.Net_nodes[node_i].child[child_i]->parent2->label!=my_Net.Net_nodes[node_i].label){
						throw std::invalid_argument("parent2 is not correct!!!");
					}
					for (unsigned int child_i_contains_gt_node1_i=0;child_i_contains_gt_node1_i<my_Net.Net_nodes[node_i].child[child_i]->Net_node_contains_gt_node2.size();child_i_contains_gt_node1_i++){
						my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.push_back(my_Net.Net_nodes[node_i].child[child_i]->Net_node_contains_gt_node2[child_i_contains_gt_node1_i]);
					}					
				}
			}

			dout<<"In population: "<< my_Net.Net_nodes[node_i].label<<" with rank  "<<rank_i;
			
			if (rank_i==my_Net.max_rank){dout<<", at the root, everything coalesces"<<endl;}else{dout<<endl;}
			//here to choose go left or right for hybrid node.			
			if (my_Net.Net_nodes[node_i].parent2){
				dout<<"hybrid node"<<endl;
				vector <unsigned int> Net_node_contains_gt_node_dummy=my_Net.Net_nodes[node_i].Net_node_contains_gt_node1;
				my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.clear();
				double left_para=extract_hybrid_para(my_Net.Net_nodes[node_i].label);
				for (unsigned int contains_gt_node_i=0;contains_gt_node_i<Net_node_contains_gt_node_dummy.size();contains_gt_node_i++){
					if (unifRand()<left_para){
						my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.push_back(Net_node_contains_gt_node_dummy[contains_gt_node_i]);
					}
					else{
						my_Net.Net_nodes[node_i].Net_node_contains_gt_node2.push_back(Net_node_contains_gt_node_dummy[contains_gt_node_i]);
					}
				}
			}
			
			double num_lineage=double(my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.size());
						
			vector < vector <double> > lambda_bk_mat;
			if (my_para_net.Net_nodes[node_i].brchlen1!=2){
				lambda_bk_mat=build_lambda_bk_mat(my_para_net.Net_nodes[node_i].brchlen1,num_lineage);			
			}
			
			double remaining_length=my_Net.Net_nodes[node_i].brchlen1;
			dout<<num_lineage<<" lineages entering population "<< my_Net.Net_nodes[node_i].label <<" with remaining branch length of "<<remaining_length  <<endl;

			int nc=0;
			double length=0;
			if (num_lineage>1 && nc <2){
				if (my_para_net.Net_nodes[node_i].brchlen1!=2){
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
			unsigned int gt_child_node_index;
			while (((length < remaining_length) && (num_lineage>1) )  || ( (rank_i==my_Net.max_rank) && (num_lineage>1 )) ){
					
					//dout<<endl<<endl;
					//dout<<"(1)num_lineage "<<num_lineage<<endl;
					//dout<<"(1)remaining length is   "<<remaining_length<<endl;
					//dout<<"(1)proposed exp length is   "<<length<<endl;		
				if (nc>1){
					dout<<"  "<<nc<<" lineages coalesce at time "<<length<<endl;
					for (int nc_i=0;nc_i<nc;nc_i++){
						lineage_index= rand() % my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.size();
						gt_child_node_index=my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[lineage_index];
						gt_nodes_ptr[gt_child_node_index]->brchlen1=gt_nodes_ptr[gt_child_node_index]->brchlen1+length; // b_alpha <- b_alpha + l

						add_node(gt_nodes_ptr[remaining_gt_node[0]],gt_nodes_ptr[gt_child_node_index]);
						if (sim_num_gener_bool){
							num_gener_gt_nodes_ptr[gt_child_node_index]->brchlen1=num_gener_gt_nodes_ptr[gt_child_node_index]->brchlen1+(length * my_pop_net.Net_nodes[node_i].brchlen1) ; // b_alpha <- b_alpha + l*pop_size
							add_node(num_gener_gt_nodes_ptr[remaining_gt_node[0]],num_gener_gt_nodes_ptr[gt_child_node_index]);
						}
						my_gt_coal_unit.descndnt[remaining_gt_node[0]]=my_gt_coal_unit.descndnt[remaining_gt_node[0]]+my_gt_coal_unit.descndnt[gt_child_node_index];
						if (sim_num_gener_bool){
							my_gt_num_gener.descndnt[remaining_gt_node[0]]=my_gt_num_gener.descndnt[remaining_gt_node[0]]+my_gt_num_gener.descndnt[gt_child_node_index];
						}
						my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.erase(my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.begin()+lineage_index); // X'\{alpha}
					}
				}
				
				my_gt_coal_unit.Net_nodes[remaining_gt_node[0]].num_descndnt=my_gt_coal_unit.descndnt[remaining_gt_node[0]].sum();
				if (sim_num_gener_bool){
					my_gt_num_gener.Net_nodes[remaining_gt_node[0]].num_descndnt=my_gt_num_gener.descndnt[remaining_gt_node[0]].sum();
				}
				my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.push_back(remaining_gt_node[0]);// introducing a new lineage gamma


				gt_nodes_ptr[remaining_gt_node[0]]->absolute_time=gt_nodes_ptr[remaining_gt_node[0]]->child[0]->absolute_time+gt_nodes_ptr[remaining_gt_node[0]]->child[0]->brchlen1; // a_gamma <- a_alpha + b_alpha
				for (unsigned int update_brch_i=0;update_brch_i<my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.size()-1;update_brch_i++){
					gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[update_brch_i]]->brchlen1=gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[update_brch_i]]->brchlen1+length;		
					if (sim_num_gener_bool){
						num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[update_brch_i]]->brchlen1= length*my_pop_net.Net_nodes[node_i].brchlen1+num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[update_brch_i]]->brchlen1;
					}
				}
				remaining_gt_node.erase(remaining_gt_node.begin());
				
				remaining_length= remaining_length - length;
					//dout<<"remaining_length "<<remaining_length<<endl;
				num_lineage=double(my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.size());

					//dout<<"checking num_lineage "<<num_lineage<<endl;

				nc=0;
				if (num_lineage>1 && nc <2){
					if (my_para_net.Net_nodes[node_i].brchlen1!=2){
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
				dout<<num_lineage<<" live lineages in population "<< my_Net.Net_nodes[node_i].label <<" with remaining branch length of "<<remaining_length  <<endl;
			
			}
			
				
			
			if (my_Net.Net_nodes[node_i].parent2){			
				dout<<"hybrid node parent1 finished, in parent 2 now"<<endl;
				double num_lineage=(my_Net.Net_nodes[node_i].Net_node_contains_gt_node2.size());
				lambda_bk_mat.clear();
				if (my_para_net.Net_nodes[node_i].brchlen2!=2){
					lambda_bk_mat=build_lambda_bk_mat(my_para_net.Net_nodes[node_i].brchlen2,num_lineage);
				}
											
				double remaining_length=my_Net.Net_nodes[node_i].brchlen2;
				dout<<num_lineage<<" lineages entering population "<< my_Net.Net_nodes[node_i].label <<" with remaining branch length of "<<remaining_length  <<endl;

				double length=0;		
				int nc=0;
				if (num_lineage>1 && nc <2){
					if (my_para_net.Net_nodes[node_i].brchlen2!=2){
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
				unsigned int gt_child_node_index;
				while ( (length < remaining_length && num_lineage>1.0 )  ){
						//dout<<endl<<endl;
						//dout<<"(1)num_lineage "<<num_lineage<<endl;
						//dout<<"(1)remaining length is   "<<remaining_length<<endl;
						//dout<<"(1)proposed exp length is   "<<length<<endl;		

					if (nc>1){
						dout<<"  "<<nc<<" lineages coalesce at time "<<length<<endl;

						for (int nc_i=0;nc_i<nc;nc_i++){
							lineage_index= rand() % my_Net.Net_nodes[node_i].Net_node_contains_gt_node2.size();
							gt_child_node_index=my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[lineage_index];
							gt_nodes_ptr[gt_child_node_index]->brchlen1=gt_nodes_ptr[gt_child_node_index]->brchlen1+length; // b_alpha <- b_alpha + l
							add_node(gt_nodes_ptr[remaining_gt_node[0]],gt_nodes_ptr[gt_child_node_index]);
							if (sim_num_gener_bool){
								num_gener_gt_nodes_ptr[gt_child_node_index]->brchlen1=num_gener_gt_nodes_ptr[gt_child_node_index]->brchlen1+length*my_pop_net.Net_nodes[node_i].brchlen2; // b_alpha <- b_alpha + l
								add_node(num_gener_gt_nodes_ptr[remaining_gt_node[0]],num_gener_gt_nodes_ptr[gt_child_node_index]);
							}
							my_gt_coal_unit.descndnt[remaining_gt_node[0]]=my_gt_coal_unit.descndnt[remaining_gt_node[0]]+my_gt_coal_unit.descndnt[gt_child_node_index];
							if (sim_num_gener_bool){
								my_gt_num_gener.descndnt[remaining_gt_node[0]]=my_gt_num_gener.descndnt[remaining_gt_node[0]]+my_gt_num_gener.descndnt[gt_child_node_index];
							}
							my_Net.Net_nodes[node_i].Net_node_contains_gt_node2.erase(my_Net.Net_nodes[node_i].Net_node_contains_gt_node2.begin()+lineage_index); // X'\{alpha}
							
						}
					}

					my_gt_coal_unit.Net_nodes[remaining_gt_node[0]].num_descndnt=my_gt_coal_unit.descndnt[remaining_gt_node[0]].sum();
					if (sim_num_gener_bool){
						my_gt_num_gener.Net_nodes[remaining_gt_node[0]].num_descndnt=my_gt_num_gener.descndnt[remaining_gt_node[0]].sum();
					}
					my_Net.Net_nodes[node_i].Net_node_contains_gt_node2.push_back(remaining_gt_node[0]);

					gt_nodes_ptr[remaining_gt_node[0]]->absolute_time=gt_nodes_ptr[remaining_gt_node[0]]->child[0]->absolute_time+gt_nodes_ptr[remaining_gt_node[0]]->child[0]->brchlen1;
					for (unsigned int update_brch_i=0;update_brch_i<my_Net.Net_nodes[node_i].Net_node_contains_gt_node2.size()-1;update_brch_i++){
						gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[update_brch_i]]->brchlen1=gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[update_brch_i]]->brchlen1+length;		
						if (sim_num_gener_bool){
							num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[update_brch_i]]->brchlen1= length*my_pop_net.Net_nodes[node_i].brchlen2 + num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[update_brch_i]]->brchlen1;
						}
					}
					remaining_gt_node.erase(remaining_gt_node.begin());
					remaining_length= remaining_length - length;
					num_lineage=double(my_Net.Net_nodes[node_i].Net_node_contains_gt_node2.size());
						//dout<<"checking num_lineage "<<num_lineage<<endl;
					
					nc=0;
					if (num_lineage>1 && nc <2){
						if (my_para_net.Net_nodes[node_i].brchlen2!=2){
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
					dout<<num_lineage<<" live lineages in population "<< my_Net.Net_nodes[node_i].label <<" with remaining branch length of "<<remaining_length  <<endl;

				}
			}
					
			dout<<"*************************before adjusting***************"<<endl;
			my_gt_coal_unit.print_all_node_dout()	;
			
			if (rank_i<my_Net.max_rank){
				if (my_Net.Net_nodes[node_i].parent2){
					for (unsigned int num_lineage_i=0;num_lineage_i<my_Net.Net_nodes[node_i].Net_node_contains_gt_node2.size();num_lineage_i++){
					//	gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node[num_lineage_i]]->brchlen1=gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node[num_lineage_i]]->brchlen1+my_Net.Net_nodes[node_i].parent1->absolute_time-gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node[num_lineage_i]]->absolute_time;
					//	gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node[num_lineage_i]]->absolute_time=my_Net.Net_nodes[node_i].parent1->absolute_time;	
						if (sim_num_gener_bool){
							if ((gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[num_lineage_i]]->absolute_time)>  my_Net.Net_nodes[node_i].absolute_time){
							  num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[num_lineage_i]]->brchlen1= (my_Net.Net_nodes[node_i].parent2->absolute_time - gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[num_lineage_i]]->absolute_time)*my_pop_net.Net_nodes[node_i].brchlen2;//+num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1;						
							}
							else{
							  num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[num_lineage_i]]->brchlen1= (my_Net.Net_nodes[node_i].parent2->absolute_time - gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[num_lineage_i]]->absolute_time- gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[num_lineage_i]]->brchlen1)*my_pop_net.Net_nodes[node_i].brchlen2+num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[num_lineage_i]]->brchlen1;						
							}
						}
						gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[num_lineage_i]]->brchlen1=my_Net.Net_nodes[node_i].parent2->absolute_time -gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node2[num_lineage_i]]->absolute_time;
					}					
				}
				//else{
				// through parent 1, it should always be checked!!!	
				
				for (unsigned int num_lineage_i=0;num_lineage_i<my_Net.Net_nodes[node_i].Net_node_contains_gt_node1.size();num_lineage_i++){
					if (sim_num_gener_bool){
						if ((gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[num_lineage_i]]->absolute_time)>  my_Net.Net_nodes[node_i].absolute_time){
							num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1= (my_Net.Net_nodes[node_i].parent1->absolute_time - gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[num_lineage_i]]->absolute_time)*my_pop_net.Net_nodes[node_i].brchlen1;//+num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1;						
						}
						else{
							num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1= (my_Net.Net_nodes[node_i].parent1->absolute_time - gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[num_lineage_i]]->absolute_time- gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1)*my_pop_net.Net_nodes[node_i].brchlen1+num_gener_gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1;						
						}
					}
					gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[num_lineage_i]]->brchlen1=my_Net.Net_nodes[node_i].parent1->absolute_time -gt_nodes_ptr[my_Net.Net_nodes[node_i].Net_node_contains_gt_node1[num_lineage_i]]->absolute_time;
				}
				//}
				
					dout<<"************************* after adjusting***************"<<endl;
					my_gt_coal_unit.print_all_node_dout()	;

			}
			remaining_sp_node.erase(remaining_sp_node.begin()+remaining_sp_node_i);
				//dout<<rank_i<<" "<<my_Net.Net_nodes[node_i].label<<" "<<my_Net.Net_nodes[node_i].path_time.size()<<endl;
		}
		else{
			remaining_sp_node_i++;
		}

		if (remaining_sp_node_i==remaining_sp_node.size()-1){
			rank_i++;
			remaining_sp_node_i=0;
		}
	}
	
	for (unsigned int node_i=0;node_i<gt_nodes_ptr.size();){
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
	
	
	for (unsigned int node_i=0;node_i<my_gt_coal_unit.Net_nodes.size();node_i++){
		if (my_gt_coal_unit.descndnt[node_i].sum()==0){
			break;
		}
		if (!my_gt_coal_unit.Net_nodes[node_i].tip_bool){
			for (unsigned int tax_name_i=0;tax_name_i<tax_name.size();tax_name_i++){
				if (my_gt_coal_unit.descndnt[node_i][tax_name_i] > 0 ){
					if (my_gt_coal_unit.Net_nodes[node_i].clade.size()==0){
						my_gt_coal_unit.Net_nodes[node_i].clade=tax_name[tax_name_i];
					}
					else{
						my_gt_coal_unit.Net_nodes[node_i].clade=my_gt_coal_unit.Net_nodes[node_i].clade+tax_name[tax_name_i];
					}
					my_gt_coal_unit.Net_nodes[node_i].clade.push_back('&');
				}
			}
			my_gt_coal_unit.Net_nodes[node_i].clade.erase(my_gt_coal_unit.Net_nodes[node_i].clade.size()-1,1);
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
	
	
	if (my_action.sim_mut_unit_bool){
		build_gt_string_mut_unit(mutation_rate);
	}	
	
	
	if (my_action.sim_num_mut_bool || my_action.Si_num_bool){
		Net mt_tree(gt_string_gener_num);
		vector <double> brch_total;
		vector <Node*> mt_nodes_ptr;
		total_brchlen=0;
		for (unsigned int gt_node_i=0;gt_node_i<mt_tree.Net_nodes.size();gt_node_i++){
			Node* new_node_ptr=NULL;
			mt_nodes_ptr.push_back(new_node_ptr);
			total_brchlen=total_brchlen+mt_tree.Net_nodes[gt_node_i].brchlen1;
			mt_nodes_ptr[gt_node_i]=&mt_tree.Net_nodes[gt_node_i];
			mt_nodes_ptr[gt_node_i]->brchlen1=0;
			brch_total.push_back(total_brchlen);
		}	
		
		brch_total.back()=0;
		//cout<<total_brchlen<<endl;
		//double theta=1.0;
		int total_mut = poisson_rand_var(mutation_rate*total_brchlen);
		//cout<<total_mut<<endl;
		for (int mut_i=0;mut_i<total_mut;mut_i++){
			unsigned int brch_index=0;
			double u=unifRand()*total_brchlen;
			while (u>brch_total[brch_index]){
				 brch_index =brch_index + 1;	
			}
			mt_nodes_ptr[brch_index]->brchlen1 = mt_nodes_ptr[brch_index]->brchlen1+1;
		}
		rewrite_node_content(mt_nodes_ptr);
		gt_string_mut_num=mt_nodes_ptr.back()->node_content+mt_nodes_ptr.back()->label+";";
		gt_string_mut_num=remove_interior_label(gt_string_mut_num);
		
		if (my_action.Si_num_bool){	
			Si_num_out_table(mt_tree,total_mut);
		}
	}

	if (tax_name.size()==2 && my_action.mono_bool){
		compute_monophyly_vec(my_gt_coal_unit,sample_size);
	}
	dout<<"	end of sim_one_gt::sim_one_gt(sim::param sim_param,action_board my_action)"<<endl;
}


void sim_one_gt::Si_num_out_table(Net mt_tree,int total_mut){
	vector <int> S_i(mt_tree.tip_name.size()-1,0);
	for (unsigned int sii=0;sii<S_i.size();sii++){
		for (unsigned int node_i=0;node_i<mt_tree.Net_nodes.size();node_i++){
			int sii_pluse_1=sii+1;
			if ((mt_tree.Net_nodes[node_i].brchlen1>0) && (sii_pluse_1==mt_tree.descndnt2[node_i].sum() )){
				//if ((sii+1)==mt_tree.descndnt2[node_i].sum()){
					S_i[sii]=S_i[sii]+mt_tree.Net_nodes[node_i].brchlen1;
				//}
			}
		}
	}
		
	ofstream out_table_file;
	out_table_file.open ("out_table", ios::out | ios::app | ios::binary); 
	out_table_file <<setw(12)<< total_brchlen<<setw(14)<<mt_tree.Net_nodes.back().absolute_time<<setw(7)<<total_mut<<"  ";
	for (unsigned int sii=0;sii<S_i.size();sii++){
		out_table_file<<setw(4)<<S_i[sii]<<" ";
	}
	out_table_file<<endl;
	out_table_file.close();
	
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
				if (isnan(lambda_bk_mat_b_k)){
					throw std::domain_error("Function \"build_lambda_bk_mat\" returns NaN");
					//dout<<"log(n_choose_k(b_i,k_i)) " <<log(n_choose_k(b_i,k_i))<<endl;
					//dout<<"b="<<b_i<<"  k="<<k_i<<endl;
					//dout<<"log(boost::math::binomial_coefficient(b_i,k_i)) " <<log(boost::math::binomial_coefficient<double>(unsigned(b_i),unsigned(k_i)))<<endl;
					//dout<<"log(pow(para,k_i)) "<< log(pow(para,k_i))<<endl;
					//dout<<"log(pow(1-para,b_i-k_i)) "<<log(pow(1-para,b_i-k_i))<<endl;
				}
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
	for (unsigned int gt_node_i=0;gt_node_i<gt_mut_unit.Net_nodes.size();gt_node_i++){
		Node* new_node_ptr=NULL;
		gt_mut_unit_nodes_ptr.push_back(new_node_ptr);
		gt_mut_unit_nodes_ptr[gt_node_i]=&gt_mut_unit.Net_nodes[gt_node_i];
		gt_mut_unit_nodes_ptr[gt_node_i]->brchlen1=gt_mut_unit_nodes_ptr[gt_node_i]->brchlen1 * mutation_rate;
	}
	rewrite_node_content(gt_mut_unit_nodes_ptr);
	gt_string_mut_unit=gt_mut_unit_nodes_ptr.back()->node_content+gt_mut_unit_nodes_ptr.back()->label+";";
	gt_string_mut_unit=remove_interior_label(gt_string_mut_unit);
}


void sim_one_gt::compute_monophyly_vec(Net my_gt_coal_unit,vector < int > sample_size){
	vector <double> monophyly_initial(6,0);
	monophyly=monophyly_initial;
	for (unsigned int tax_i=0;tax_i<2;tax_i++){
		for (unsigned int node_i=0;node_i<my_gt_coal_unit.Net_nodes.size();node_i++){
			if (my_gt_coal_unit.descndnt[node_i].sum()==0){
				break;
			}
			if ((my_gt_coal_unit.descndnt[node_i][tax_i]==sample_size[tax_i]) &&  ( my_gt_coal_unit.descndnt[node_i][tax_i]==my_gt_coal_unit.descndnt[node_i].sum() ) ){
				monophyly[tax_i]=1;
				break;
			}
		}
	}
	if (monophyly[0]==1){
		if (monophyly[1]==1){
			monophyly[2]=1;
		}
		else{
			monophyly[3]=1;
		}	
	}
	else{
		if (monophyly[1]==1){
			monophyly[4]=1;
		}
		else{
			monophyly[5]=1;
		}
	}

	for (size_t mono_i=0;mono_i<6;mono_i++){
		dout<<monophyly[mono_i]<<"  ";
	}
	dout<<endl;
}



/*! \brief  sim_n_gt constructor */
//sim_n_gt::sim_n_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string, vector < int > sample_size,double mutation_rate,int num_sim_gt,bool sim_mut_unit_bool, bool sim_num_gener_bool,bool sim_num_mut_bool,bool mono_bool){
//sim_n_gt::sim_n_gt(string sp_string_coal_unit, string sp_string_pop_size, string para_string, vector < int > sample_size,double mutation_rate,int num_sim_gt,action_board my_action){
sim_n_gt::sim_n_gt(sim::param sim_param,action_board my_action){
	dout<<" Start simulating "<<sim_param.num_sim_gt <<" gene trees -- begin of sim_n_gt::sim_n_gt(sim::param sim_param,action_board my_action)"<<endl;
	
	string gene_tree_file=my_action.gene_tree_file;
	string sp_string_coal_unit=sim_param.sp_string_coal_unit;
	string sp_string_pop_size=sim_param.sp_string_pop_size;
	string para_string=sim_param.para_string;
	vector < int > sample_size=sim_param.sample_size;
	double mutation_rate=sim_param.mutation_rate;
	int num_sim_gt=sim_param.num_sim_gt;
	
	ofstream sim_gt_file_coal_unit;
	string gene_tree_file_coal_unit=gene_tree_file+"_coal_unit";
	remove(gene_tree_file_coal_unit.c_str());
	sim_gt_file_coal_unit.open (gene_tree_file_coal_unit.c_str(), ios::out | ios::app | ios::binary); 
			
	ofstream sim_gt_file_mut_unit;
	ofstream sim_gt_file_num_gener;
	ofstream sim_gt_file_num_mut;
	
	string gene_tree_file_mut_unit=gene_tree_file+"_mut_unit";
	remove(gene_tree_file_mut_unit.c_str());
	string gene_tree_file_num_gener=gene_tree_file+"_num_gener";
	remove(gene_tree_file_num_gener.c_str());
	string gene_tree_file_num_mut=gene_tree_file+"_num_mut";
	remove(gene_tree_file_num_mut.c_str());
	
	if (my_action.sim_mut_unit_bool){
		sim_gt_file_mut_unit.open (gene_tree_file_mut_unit.c_str(), ios::out | ios::app | ios::binary); 
	}
	if (my_action.sim_num_gener_bool){
		sim_gt_file_num_gener.open (gene_tree_file_num_gener.c_str(), ios::out | ios::app | ios::binary); 
	}
	if (my_action.sim_num_mut_bool){
		sim_gt_file_num_mut.open (gene_tree_file_num_mut.c_str(), ios::out | ios::app | ios::binary); 
	}

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
		gt_string_coal_unit_s.push_back(sim_gt_string.gt_string_coal_unit);
		if (my_action.sim_num_mut_bool){
			gt_string_mut_num_s.push_back(sim_gt_string.gt_string_mut_num);
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
		
		//total_brchlen.push_back(sim_gt_string.total_brchlen);
		dout<<sim_gt_string.gt_string_coal_unit<<endl;
		sim_gt_file_coal_unit<<sim_gt_string.gt_string_coal_unit <<"\n";
		
		if (my_action.sim_mut_unit_bool){
			sim_gt_file_mut_unit<<sim_gt_string.gt_string_mut_unit <<"\n";
		}
		if (my_action.sim_num_gener_bool){
			sim_gt_file_num_gener<<sim_gt_string.gt_string_gener_num <<"\n";
		}
		if (my_action.sim_num_mut_bool){
			sim_gt_file_num_mut<<sim_gt_string.gt_string_mut_num <<"\n";
		}
	}
	
	if (my_action.mono_bool){
		for (unsigned int mono_i=0;mono_i<monophyly.size();mono_i++){
			monophyly[mono_i]=monophyly[mono_i]/num_sim_gt;
		}
	}
	

	sim_gt_file_coal_unit.close();
	
	if (my_action.sim_mut_unit_bool){
		sim_gt_file_mut_unit.close();
	}
	if (my_action.sim_num_gener_bool){
		sim_gt_file_num_gener.close();
	}
	if (my_action.sim_num_mut_bool){
		sim_gt_file_num_mut.close();
	}
	dout<<"end of sim_n_gt::sim_n_gt(sim::param sim_param,action_board my_action)"<<endl;

}






/*! \brief Convert the network branch lengths from number of generations to coalescent unit, by dividing the population size 
 * \return string Network in extended newick form, branch lengths are in coalescent unit*/
string write_sp_string_in_coal_unit(string sp_num_gener_string /*! Network in extended newick form, branch lengths are the number of generations*/,
string pop_size_string /*! Network in extended newick form, branch lengths are the population sizes*/){
	Net sp_num_gener_net(sp_num_gener_string);
	Net pop_size_net(pop_size_string);
	
	vector <Node*> sp_num_gener_net_node_ptr;
	for (unsigned int node_i=0;node_i<sp_num_gener_net.Net_nodes.size();node_i++){
		Node* new_node_ptr=NULL;
        sp_num_gener_net_node_ptr.push_back(new_node_ptr);
        sp_num_gener_net_node_ptr[node_i]=&sp_num_gener_net.Net_nodes[node_i];
		if (node_i<sp_num_gener_net.Net_nodes.size()-1){
			sp_num_gener_net_node_ptr[node_i]->brchlen1=(sp_num_gener_net_node_ptr[node_i]->brchlen1) / pop_size_net.Net_nodes[node_i].brchlen1;
			if (pop_size_net.Net_nodes[node_i].hybrid){
				sp_num_gener_net_node_ptr[node_i]->brchlen2=(sp_num_gener_net_node_ptr[node_i]->brchlen2) / pop_size_net.Net_nodes[node_i].brchlen2;
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
	for (unsigned int node_i=0;node_i<para_net_check.Net_nodes.size();node_i++){
		Node* new_node_ptr=NULL;
		pop_size_node_ptr.push_back(new_node_ptr);
		pop_size_node_ptr[node_i]=&pop_size_check.Net_nodes[node_i];
		if (para_net_check.Net_nodes[node_i].brchlen1<2 && para_net_check.Net_nodes[node_i].brchlen1>1){ // rescale the number of generations for alpha
			pop_size_node_ptr[node_i]->brchlen1 = pow(pop_size_node_ptr[node_i]->brchlen1,para_net_check.Net_nodes[node_i].brchlen1-1);
		}
		if (para_net_check.Net_nodes[node_i].hybrid){
			if (para_net_check.Net_nodes[node_i].brchlen2<2 && para_net_check.Net_nodes[node_i].brchlen2>1){
				pop_size_node_ptr[node_i]->brchlen2=pow(pop_size_node_ptr[node_i]->brchlen2,para_net_check.Net_nodes[node_i].brchlen2-1);
			}
		}
	}
	rewrite_node_content(pop_size_node_ptr);
	string pop_size_string_return = construct_adding_new_Net_str(pop_size_check);
	return pop_size_string_return;
}




/*! \brief Printing out the header for out_table*/
void outtable_header(int total_lineage){
	ofstream out_table_file;
	out_table_file.open ("out_table", ios::out | ios::app | ios::binary); 
	out_table_file <<"t_total       t_MRCA       S_total  ";
	for (int sii=0;sii<total_lineage-1 ;sii++){
		out_table_file<<"S_"<<sii+1<<"  ";
	}
	out_table_file<<endl;
	out_table_file.close();
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
	for (unsigned int kmerge=0;kmerge<nc_X.size();kmerge++){
		nc_X[kmerge]= -log( 1-unifRand() )/ lambda_bk_mat[num_lineage-2][kmerge];
		//cout<<nc_X[kmerge]<<endl;
		//if (isnan(nc_X[kmerge])){cout<<lambda_bk_mat[num_lineage-2][kmerge]<<endl;}
	}
	return nc_X;
}

//use heap structure for this!
int update_nc(valarray <double> nc_X){	
	for (int kmerge=0;kmerge<int(nc_X.size());kmerge++){
		if (nc_X[kmerge]== nc_X.min()){
			//nc=kmerge+2;
			return kmerge+2;
			//break;
		}
	}
}

double kingman_bl(double num_lineage){
	return -log( 1-unifRand() ) / num_lineage / (num_lineage-1) * 2.0;
}
