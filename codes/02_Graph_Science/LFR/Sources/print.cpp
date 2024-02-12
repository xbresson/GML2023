#if !defined(PRINT_INCLUDED)
#define PRINT_INCLUDED	
	


int cherr() {
	
	cerr<<"the check failed"<<endl;
	int e;
	cin>>e;
	
}


int cherr(double a) {
	
	cerr<<"the check failed because of "<<a<<endl;
	int e;
	cin>>e;
	
}



template <typename uno, typename due>
void prints(pair <uno, due> &sq,  ostream &out) {

	out<<sq.first<<"\t"<<sq.second<<endl;
	
}


template <typename uno, typename due>
void prints(pair <uno, due> &sq) {

	cout<<sq.first<<"\t"<<sq.second<<endl;
	
}


template <typename uno, typename due>
void prints(map <uno, due> &sq,  ostream &out) {

	typename map <uno, due>::iterator it = sq.begin(); 
	while(it != sq.end()) { 
		out<<it->first<<"\t"<<it->second<<endl;
		it++; 
	} 

	out<<endl;
	
}

template <typename uno, typename due>
void prints(multimap <uno, due> &sq,  ostream &out) {

	typename map <uno, due>::iterator it = sq.begin(); 
	while(it != sq.end()) { 
		out<<it->first<<"\t"<<it->second<<endl;
		it++; 
	} 

	out<<endl;
	
}



template <typename Seq>
void prints(Seq &sq, ostream &out) {

	typename Seq::iterator it = sq.begin(); 
	while(it != sq.end())
		out<<*(it++)<<"\t";
		

	out<<endl;
	
}

template <typename type_>
void prints(type_ *a, int b) {
	
	for (int i=0; i<b; i++)
		cout<<a[i]<<" ";
	cout<<endl;


}


template<typename T, template<typename> class C>
void printm(C<T>& c, ostream &out) { 
	
	typename C<T>::iterator it = c.begin(); 
	while(it != c.end()) { 
		prints(*it, out);
		it++; 
	} 

	out<<endl;
} 
 
	


template <typename uno, typename due>
void prints(map <uno, due> &sq) {

	typename map <uno, due>::iterator it = sq.begin(); 
	while(it != sq.end()) { 
		cout<<it->first<<"\t"<<it->second<<endl;
		it++; 
	} 

	cout<<endl;
	
}

template <typename uno, typename due>
void prints(multimap <uno, due> &sq) {

	typename map <uno, due>::iterator it = sq.begin(); 
	while(it != sq.end()) { 
		cout<<it->first<<"\t"<<it->second<<endl;
		it++; 
	} 

	cout<<endl;
	
}



template <typename Seq>
void prints(Seq &sq) {

	typename Seq::iterator it = sq.begin(); 
	while(it != sq.end())
		cout<<*(it++)<<"\t";
		

	cout<<endl;
	
}


template <typename type>
void prints(const deque<type> & sq) {
	
	for(int i=0; i<sq.size(); i++)
		cout<<sq[i]<<"\t";
	cout<<endl;
	
		
}





template<typename T, template<typename> class C>
void printm(C<T>& c) { 
	
	typename C<T>::iterator it = c.begin(); 
	while(it != c.end()) { 
		prints(*it);
		it++; 
	} 

	cout<<endl;
} 







#endif
