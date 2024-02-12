


#if !defined(HISTOGRAMS_INCLUDED)
#define HISTOGRAMS_INCLUDED








template <typename type>
int log_histogram(deque<type> &c, ostream & out, int number_of_bins) {		// c is the set od data, min is the lower bound, max is the upper one
	
	
	
	deque <type> d;
	for(int i=0; i<c.size(); i++) if (c[i]>0)
		d.push_back(c[i]);
	
	c.clear();
	c=d;
	
	double min=double(c[0]);
	double max=double(c[0]);
	
	for (int i=0; i<c.size(); i++) {
		
		if (min>double(c[i]))
			min=double(c[i]);
		
		if (max<double(c[i]))
			max=double(c[i]);
		
	}
	
	
	
	
	deque <int> hist;
	deque <double> hist2;
	deque <double> bins;
	double step=log(min);
	if (max==min)
		max++;
	
	double bin=(log(max)-log(min))/number_of_bins;		// bin width
	
		

	while (step<=log(max)+2*bin) {
		
		
		bins.push_back(exp(step));
		hist.push_back(0);			
		hist2.push_back(0);			
		step+=bin;
	}
	
	
	for (int i=0; i<c.size(); i++) {
		
		
		int index=bins.size()-1;
		for (int j=0; j<bins.size()-1; j++) if(	(fabs(double(c[i])-bins[j])<1e-7) || (	double(c[i])>bins[j]	&&	double(c[i])<bins[j+1]	)	) { 
		// this could be done in a more efficient way
			
			index=j;
			break;
		
		}
		
		//cout<<hist[index]<<" "<<index<<endl;
		
				
		hist[index]++;
		hist2[index]+=double(c[i]);
		
	}
	
	
	
	
	for (int i=0; i<hist.size()-1; i++) {
		
		double h1= bins[i];
		double h2= bins[i+1];
		
		
		double x=hist2[i]/hist[i];
		double y=double(hist[i])/(c.size()*(h2-h1));
		
		if (fabs(y)>1e-10)
			out<<x<<"\t"<<y<<endl;		
		
	
	
	}
	
	
	
	return 0;

}








template <typename type>
int histogram (vector <type> &c, ostream & out, int number_of_bins, double b1, double b2) {		

	// this should be OK
	// c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
	
	
	
	double min=double(c[0]);
	double max=double(c[0]);
	
	for (int i=0; i<c.size(); i++) {
		
		if (min>double(c[i]))
			min=double(c[i]);
		
		if (max<double(c[i]))
			max=double(c[i]);
		
	}
	
	
	
	min-=1e-6;
	max+=1e-6;
	
	
	
	if (b1!=b2) {
		
		min=b1;
		max=b2;
	
	}
		
	if (max==min)
		max+=1e-3;
	
	
	
	deque <int> hist;
	deque <double> hist2;
		
	double step=min;
	double bin=(max-min)/number_of_bins;		// bin width

	while (step<=max+2*bin) {
	
		hist.push_back(0);			
		hist2.push_back(0);			
		step+=bin;
	}
	

	
		
	
	
	for (int i=0; i<c.size(); i++) {
		
		
		
		double data=double(c[i]);
		
		if (data>min && data<=max) {
			
			int index=int((data-min)/bin);		
			
				
			hist[index]++;
			hist2[index]+=double(c[i]);
		
		}
		
	}
	
	for (int i=0; i<hist.size()-1; i++) {
		
		
		
				
		double x=hist2[i]/hist[i];
		double y=double(hist[i])/(c.size()*bin);
		
		if (fabs(y)>1e-10)
			out<<x<<"\t"<<y<<endl;
		
	
	}
	
	
	
			
	return 0;

}





template <typename type>
int not_norm_histogram_correlated (deque<type> &c, deque<type> &d, ostream & out, int number_of_bins, double b1, double b2) {		
	
	
	// c is the x axis, d the y, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
	
	
	
	double min=double(c[0]);
	double max=double(c[0]);
	
	for (int i=0; i<c.size(); i++) {
		
		if (min>double(c[i]))
			min=double(c[i]);
		
		if (max<double(c[i]))
			max=double(c[i]);
		
	}
	
	
	
	min-=1e-6;
	max+=1e-6;
	
	
	
	if (b1!=b2) {
		
		min=b1;
		max=b2;
	
	}
		
	if (max==min)
		max+=1e-3;
	
	
	
	deque <int> hist;			// frequency in the bin
	deque <double> hist_x;		// x sum in the bin
	deque <double> hist_y;		// y sum in the bin
		
	double step=min;
	double bin=(max-min)/number_of_bins;		// bin width

	while (step<=max+2*bin) {
	
		hist.push_back(0);			
		hist_x.push_back(0);			
		hist_y.push_back(0);			
		step+=bin;
	}
	

	
		
	
	
	for (int i=0; i<c.size(); i++) {
		
		
		
		double data=double(c[i]);
		
		if (data>min && data<=max) {
			
			int index=int((data-min)/bin);		
			
				
			hist[index]++;
			hist_x[index]+=double(c[i]);
			hist_y[index]+=double(d[i]);
		
		}
		
	}
	
	for (int i=0; i<hist.size()-1; i++) {
		
		
		
				
		double x=hist_x[i]/hist[i];
		double y=hist_y[i]/hist[i];;
		
		if (fabs(y)>1e-10)
			out<<x<<"\t"<<y<<endl;
		
	
	}
	
	
	
			
	return 0;

}




template <typename type>
int histogram (deque <type> &c, ostream & out, int number_of_bins, double b1, double b2) {		

	// this should be OK
	// c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
	
	
	
	double min=double(c[0]);
	double max=double(c[0]);
	
	for (int i=0; i<c.size(); i++) {
		
		if (min>double(c[i]))
			min=double(c[i]);
		
		if (max<double(c[i]))
			max=double(c[i]);
		
	}
	
	
	
	min-=1e-6;
	max+=1e-6;
	
	
	
	if (b1!=b2) {
		
		min=b1;
		max=b2;
	
	}
		
	if (max==min)
		max+=1e-3;
	
	
	
	deque <int> hist;
	deque <double> hist2;
		
	double step=min;
	double bin=(max-min)/number_of_bins;		// bin width

	while (step<=max+2*bin) {
	
		hist.push_back(0);			
		hist2.push_back(0);			
		step+=bin;
	}
	

	
		
	
	
	for (int i=0; i<c.size(); i++) {
		
		
		
		double data=double(c[i]);
		
		if (data>min && data<=max) {
			
			int index=int((data-min)/bin);		
			
				
			hist[index]++;
			hist2[index]+=double(c[i]);
		
		}
		
	}
	
	for (int i=0; i<hist.size()-1; i++) {
		
		
		
				
		double x=hist2[i]/hist[i];
		double y=double(hist[i])/(c.size()*bin);
		
		if (fabs(y)>1e-10)
			out<<x<<"\t"<<y<<endl;
		
	
	}
	
	
	
			
	return 0;

}



template <typename type>
int not_norm_histogram (vector<type> &c, ostream & out, int number_of_bins, double b1, double b2) {		

	// this should be OK
	// c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
	
	
	
	double min=double(c[0]);
	double max=double(c[0]);
	
	for (int i=0; i<c.size(); i++) {
		
		if (min>double(c[i]))
			min=double(c[i]);
		
		if (max<double(c[i]))
			max=double(c[i]);
		
	}
	
	
	
	min-=1e-6;
	max+=1e-6;
	
	
	
	if (b1!=b2) {
		
		min=b1;
		max=b2;
	
	}
		
	if (max==min)
		max+=1e-3;
	
	
	
	deque <int> hist;
	deque <double> hist2;
		
	double step=min;
	double bin=(max-min)/number_of_bins;		// bin width

	while (step<=max+2*bin) {
	
		hist.push_back(0);			
		hist2.push_back(0);			
		step+=bin;
	}
	

	
		
	
	
	for (int i=0; i<c.size(); i++) {
		
		
		
		double data=double(c[i]);
		
		if (data>min && data<=max) {
			
			int index=int((data-min)/bin);		
			
				
			hist[index]++;
			hist2[index]+=double(c[i]);
		
		}
		
	}
	
	
	for (int i=0; i<hist.size()-1; i++) {
		
		
		
				
		double x=hist2[i]/hist[i];
		double y=double(hist[i])/(c.size());
		
		if (fabs(y)>1e-10)
			out<<x<<"\t"<<y<<endl;
		
	
	}
	
	
	
			
	return 0;

}



template <typename type>
int not_norm_histogram (deque<type> &c, ostream & out, int number_of_bins, double b1, double b2) {		

	// this should be OK
	// c is the set of data, b1 is the lower bound, b2 is the upper one (if they are equal, default limits are used)
	
	
	
	double min=double(c[0]);
	double max=double(c[0]);
	
	for (int i=0; i<c.size(); i++) {
		
		if (min>double(c[i]))
			min=double(c[i]);
		
		if (max<double(c[i]))
			max=double(c[i]);
		
	}
	
	
	
	min-=1e-6;
	max+=1e-6;
	
	
	
	if (b1!=b2) {
		
		min=b1;
		max=b2;
	
	}
		
	if (max==min)
		max+=1e-3;
	
	
	
	deque <int> hist;
	deque <double> hist2;
		
	double step=min;
	double bin=(max-min)/number_of_bins;		// bin width

	while (step<=max+2*bin) {
	
		hist.push_back(0);			
		hist2.push_back(0);			
		step+=bin;
	}
	

	
		
	
	
	for (int i=0; i<c.size(); i++) {
		
		
		
		double data=double(c[i]);
		
		if (data>min && data<=max) {
			
			int index=int((data-min)/bin);		
			
				
			hist[index]++;
			hist2[index]+=double(c[i]);
		
		}
		
	}
	
	
	for (int i=0; i<hist.size()-1; i++) {
		
		
		
				
		double x=hist2[i]/hist[i];
		double y=double(hist[i])/(c.size());
		
		if (fabs(y)>1e-10)
			out<<x<<"\t"<<y<<endl;
		
	
	}
	
	
	
			
	return 0;

}

int int_histogram (vector <int> &c, ostream & out) {

	
	
	map<int, double> hist;
	
	double freq=1/double(c.size());
	
	for (int i=0; i<c.size(); i++) {
		
		map<int, double>::iterator itf=hist.find(c[i]);
		if (itf==hist.end())
			hist.insert(make_pair(c[i], 1.));
		else
			itf->second++;
	
	
	}
	
	
	for (map<int, double>::iterator it=hist.begin(); it!=hist.end(); it++)
		it->second=it->second*freq;
	
	prints(hist, out);



}

int int_histogram (deque <int> &c, ostream & out) {

	
	
	map<int, double> hist;
	
	double freq=1/double(c.size());
	
	for (int i=0; i<c.size(); i++) {
		
		map<int, double>::iterator itf=hist.find(c[i]);
		if (itf==hist.end())
			hist.insert(make_pair(c[i], 1.));
		else
			itf->second++;
	
	
	}
	
	
	for (map<int, double>::iterator it=hist.begin(); it!=hist.end(); it++)
		it->second=it->second*freq;
	
	prints(hist, out);



}




#endif


