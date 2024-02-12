#if !defined(RANDOM_INCLUDED)
#define RANDOM_INCLUDED	
	


#define R2_IM1 2147483563
#define R2_IM2 2147483399
#define R2_AM (1.0/R2_IM1)
#define R2_IMM1 (R2_IM1-1)
#define R2_IA1 40014
#define R2_IA2 40692
#define R2_IQ1 53668
#define R2_IQ2 52774
#define R2_IR1 12211
#define R2_IR2 3791
#define R2_NTAB 32
#define R2_NDIV (1+R2_IMM1/R2_NTAB)
#define R2_EPS 1.2e-7
#define R2_RNMX (1.0-R2_EPS)




double ran2(long *idum) {
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[R2_NTAB];
	double temp;

	if(*idum<=0 || !iy){
		if(-(*idum)<1) *idum=1*(*idum);
		else *idum=-(*idum);
		idum2=(*idum);
		for(j=R2_NTAB+7;j>=0;j--){
			k=(*idum)/R2_IQ1;
			*idum=R2_IA1*(*idum-k*R2_IQ1)-k*R2_IR1;
			if(*idum<0) *idum+=R2_IM1;
			if(j<R2_NTAB) iv[j]=*idum;
		}
		iy=iv[0];
	}
	k=(*idum)/R2_IQ1;
	*idum=R2_IA1*(*idum-k*R2_IQ1)-k*R2_IR1;
	if(*idum<0) *idum+=R2_IM1;
	k=(idum2)/R2_IQ2;
	idum2=R2_IA2*(idum2-k*R2_IQ2)-k*R2_IR2;
	if (idum2 < 0) idum2 += R2_IM2;
	j=iy/R2_NDIV;
	iy=iv[j]-idum2;
	iv[j]=*idum;
	if(iy<1) iy+=R2_IMM1;
	if((temp=R2_AM*iy)>R2_RNMX) return R2_RNMX;
	else return temp;
}



double ran4(bool t, long s) {
	
	double r=0;
	
	
	static long seed_=1;
	
	if(t)
		r=ran2(&seed_);
	else
		seed_=s;
	

	return r;
}


double ran4() {
	
	return ran4(true, 0);
}


void srand4(void) {
	
	long s=(long)time(NULL);
	ran4(false, s);
	
	
	
}

void srand5(int rank) {
	
	long s=(long)(rank);
	ran4(false, s);
	
}



int irand(int n) {

	return (int(ran4()*(n+1)));
	
}


void srand_file(void) {

	ifstream in("time_seed.dat");
	int seed;
	
	if (!in.is_open())
		seed=21111983;
	else
		in>>seed;
	
	if (seed < 1 || seed>R2_IM2)
		seed=1;
	
	
	srand5(seed);
	ofstream out("time_seed.dat");
	out<<seed+1<<endl;
	

}




int configuration_model(deque<set<int> > & en, deque<int> & degrees) {
	
	
	// this function is to build a network with the degree seq in degrees which is sorted (correspondence is based on the vectorial index)
	if(degrees.size()<3) {
		
		cerr<<"it seems that some communities should have only 2 nodes! This does not make much sense (in my opinion) Please change some parameters!"<<endl;
		return -1;
	
	}
	
	
	sort(degrees.begin(), degrees.end());
	

	{
		set<int> first;
		for(int i=0; i<degrees.size(); i++) 
			en.push_back(first);
	}
	
	
	
	multimap <int, int> degree_node;
	
	for(int i=0; i<degrees.size(); i++)
		degree_node.insert(degree_node.end(), make_pair(degrees[i], i));
	
	int var=0;

	while (degree_node.size() > 0) {
		
		multimap<int, int>::iterator itlast= degree_node.end();
		itlast--;
		
		multimap <int, int>::iterator itit= itlast;
		deque <multimap<int, int>::iterator> erasenda;
		
		int inserted=0;
		
		for (int i=0; i<itlast->first; i++) {
			
			if(itit!=degree_node.begin()) {
			
				itit--;
				
				
				en[itlast->second].insert(itit->second);
				en[itit->second].insert(itlast->second);
				inserted++;
				
				erasenda.push_back(itit);				
				
			}
			
			else
				break;
		
		}
		
		
		for (int i=0; i<erasenda.size(); i++) {
			
			
			if(erasenda[i]->first>1)
				degree_node.insert(make_pair(erasenda[i]->first - 1, erasenda[i]->second));
	
			degree_node.erase(erasenda[i]);
		
		}

		
		var+= itlast->first - inserted;
		degree_node.erase(itlast);
		
	}

	
	
	// this is to randomize the subgraph -------------------------------------------------------------------
	
	for(int node_a=0; node_a<degrees.size(); node_a++) for(int krm=0; krm<en[node_a].size(); krm++) {
	
					
				
		int random_mate=irand(degrees.size()-1);
		while (random_mate==node_a)
			random_mate=irand(degrees.size()-1);
				
		
		if (en[node_a].insert(random_mate).second) {
			
			deque <int> out_nodes;
			for (set<int>::iterator it_est=en[node_a].begin(); it_est!=en[node_a].end(); it_est++) if ((*it_est)!=random_mate)
				out_nodes.push_back(*it_est);
						
										
					
			int old_node=out_nodes[irand(out_nodes.size()-1)];
					
			en[node_a].erase(old_node);
			en[random_mate].insert(node_a);
			en[old_node].erase(node_a);

										
			deque <int> not_common;
			for (set<int>::iterator it_est=en[random_mate].begin(); it_est!=en[random_mate].end(); it_est++)
				if ((old_node!=(*it_est)) && (en[old_node].find(*it_est)==en[old_node].end()))
					not_common.push_back(*it_est);
					
						
			int node_h=not_common[irand(not_common.size()-1)];
			
			en[random_mate].erase(node_h);
			en[node_h].erase(random_mate);
			en[node_h].insert(old_node);
			en[old_node].insert(node_h);
			
			
		}
		
		
	}

	
	

	return 0;




}




#endif
