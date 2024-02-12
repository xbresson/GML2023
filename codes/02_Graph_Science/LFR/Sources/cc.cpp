
int common_neighbors(int a, int b, deque<set<int> > & en) {
	
	if(en[a].size()>en[b].size())
		return common_neighbors(b, a, en);
	
	int number_of_triangles=0;
	
	for (set<int>::iterator iti=en[a].begin(); iti!=en[a].end(); iti++)
		if(en[b].find(*iti)!=en[b].end())
			number_of_triangles++;

	
	
	return number_of_triangles;



}


//*
double compute_cc(deque<set<int> > & en, int i) {


		
	double number_of_triangles=0;
	for (set<int>::iterator iti=en[i].begin(); iti!=en[i].end(); iti++) {
		number_of_triangles+=common_neighbors(i, *iti, en);
		
		
	}
		
	return number_of_triangles/((en[i].size())*(en[i].size()-1.));
	
	
}



double compute_cc(deque<set<int> > & en) {


	double cc=0;
	

	for(int i=0; i<en.size(); i++) {
		
		
		double number_of_triangles=0;
		for (set<int>::iterator iti=en[i].begin(); iti!=en[i].end(); iti++) {
			number_of_triangles+=common_neighbors(i, *iti, en);
		
		}
			
		
		cc+=number_of_triangles/((en[i].size())*(en[i].size()-1.));
		
	}
	
	cc/=en.size();
	
	
	
	
	return cc;




}



double compute_tot_t(deque<set<int> > & en) {


	double number_of_triangles=0;
	

	for(int i=0; i<en.size(); i++)
		for (set<int>::iterator iti=en[i].begin(); iti!=en[i].end(); iti++)
			number_of_triangles+=common_neighbors(i, *iti, en);
	

	return number_of_triangles;



}


int choose_the_least(deque<set<int> > & en, deque<int> & A, int a, int & cn_a_o) {
	
	
	int old_node;
	shuffle_s(A);
	
	cn_a_o=en[a].size();
				
	for(int i=0; i<A.size(); i++) {
		
		int nec=common_neighbors(a, A[i], en);
		if(nec < cn_a_o) {
		
			old_node=A[i];
			cn_a_o=nec;
		}
		
		if(cn_a_o==0)
			return old_node;
	}


	return old_node;


}




int cclu(deque<set<int> > & en, const deque<deque<int> > & member_list, const deque<deque<int> > & member_matrix, double ca) {

	
		
	
	double cc0=compute_cc(en);
	cout<<"Average Clustering coefficient... "<<cc0<<" trying to reach "<<ca<<endl;
	
	
	deque<double> ccs;
	for(int i=0; i<en.size(); i++)
		ccs.push_back(compute_cc(en, i));
	
	
	
	double min_relative_inc=1e-6;
	//int number_of_triangles=compute_tot_t(en);
	
	int num_p=min(int(en.size()/10), 5);
	
	while(cc0 < ca) {
	
		
		double ccold=cc0;
		
		
		for(int y=0; y<num_p; y++) for(int Ai=0; Ai<en.size(); Ai++) {
			
			
			
			// ************************************************  rewiring
			
			while(true) {
				
				int random_node = irand(en.size()-1);
				int a=random_from_set(en[random_node]);
				
				deque<int> not_ra;
				for (set<int>::iterator it_est=en[random_node].begin(); it_est!=en[random_node].end(); it_est++) if(en[a].find(*it_est)==en[a].end() && *it_est!=a)
					not_ra.push_back(*it_est);
				
				if(not_ra.size()==0)
					break;
				
				
				int random_mate=not_ra[irand(not_ra.size()-1)];
				
				bool b1=they_are_mate(a, random_mate, member_list);
				

							
				
				deque <int> out_nodes;
				for (set<int>::iterator it_est=en[a].begin(); it_est!=en[a].end(); it_est++) if(they_are_mate(a, *it_est, member_list)==b1)
					out_nodes.push_back(*it_est);
				
				if(out_nodes.size()==0)
					break;
				
				int t1;
				int old_node = choose_the_least(en, out_nodes, a, t1);
				
				//int old_node=out_nodes[irand(out_nodes.size()-1)];
				
				deque<int> not_common;
				for (set<int>::iterator it_est=en[random_mate].begin(); it_est!=en[random_mate].end(); it_est++)
					if ((old_node!=(*it_est)) && (en[old_node].find(*it_est)==en[old_node].end())) if(they_are_mate(*it_est, random_mate, member_list)==b1  && they_are_mate(*it_est, old_node, member_list)==b1)
						not_common.push_back(*it_est);
						
				if(not_common.size()==0)
					break;
				
				
				//int node_h=not_common[irand(not_common.size()-1)];
				int t2;
				int node_h = choose_the_least(en, not_common, random_mate, t2);
				
				
				//double c1=common_neighbors(a, old_node, en) + common_neighbors(random_mate, node_h, en);
				double c1=t1 + t2;
				
				
				
				en[a].erase(old_node);
				en[a].insert(random_mate);
				
				en[old_node].erase(a);
				en[old_node].insert(node_h);
				
				en[random_mate].erase(node_h);
				en[random_mate].insert(a);
				
				
				en[node_h].erase(random_mate);
				en[node_h].insert(old_node);
				
				
				
				double c2=common_neighbors(a, random_mate, en) + common_neighbors(old_node, node_h, en);
				
			
				
				if(c1>c2) {
					
					en[a].insert(old_node);
					en[a].erase(random_mate);
				
					en[old_node].insert(a);
					en[old_node].erase(node_h);
				
					en[random_mate].insert(node_h);
					en[random_mate].erase(a);
				
				
					en[node_h].insert(random_mate);
					en[node_h].erase(old_node);
					
				
				} 
				
				
				
				break;
				
				
				
			}
			
			// ************************************************  rewiring
			
			
			
			
		}
		
		
		cc0=compute_cc(en);
		
		if(cc0-ccold < min_relative_inc * cc0) {
			
			cout<<"It seems I cannot reach the wished value. I'll stop here..."<<endl;
			break;
		
		
		}
		
		
		num_p=cast_int((ca-cc0)/ (cc0-ccold)) * num_p;
		
		if(num_p<=0)
			num_p=1;
		if(num_p>50)
			num_p=50;

		
		cout<<"Average Clustering coefficient... "<<cc0<<" trying to reach "<<ca<<"\t\t expected "<<num_p<<" more step(s) "<<endl;

	}
	

	
	
	
	return 0;



}


