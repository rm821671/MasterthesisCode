
const int Nv = 10;

template<class T, size_t N>
size_t size(T (&)[N]) { return N; }

void vectors(){
	
	TVector3 vTemp[Nv];
	TLorentzVector lvTemp[Nv];
	
	for(int i=0;i<Nv;i++){
		vTemp[i].SetXYZ(0.,0.,0.);
		lvTemp[i].SetPxPyPzE(0.,0.,0.,0.);
	}
	
	
	
	vector<float> vec;
	cout << vec.size() << endl;
	
	for(int i=0;i<5;i++){
		vec.push_back(i*i);
	}
	cout << vec.size() << endl;
	cout << vec[5] << endl;
	
}

void strings(){
	
	string tname = "abcdefghij";
	cout << tname.find("efg") << endl;
	cout << tname.length() << endl;
	cout << string::npos << endl;
	
}





void arrays(){
	
	int comb = 0;
	int j[5] = {4,5,3,7,8};
	const int N = size(j);
	
	for(int i=0,n=size(j); i<n; i++){
		
		for(int k=i+1; k<n; k++){
			comb++;
			cout << comb << " with: " << j[i] << " and " << j[k] << endl;
		}
	}
	
	
	cout << "index of max element: " << distance(j, max_element(j, j + N))<< endl;
	
	
}


void byreferencedaughter(float& a, float& b){
	
	a = a/2.;
	
	b = a+b;
	
}

void byreference(){
	
	float k = 3;
	float z = 5;
	
	byreferencedaughter(k, z);
	
	cout << k << endl;
	cout << z << endl;
	
	d = 4.5e-3;
	
	cout << d << endl;
	
}

void modulos(){
	
	int N=10;
	
	cout << N/5 << endl;
	cout << N % 5 << endl;
	cout << N % 10 << endl;
	cout << N % 3 << endl;
	
	
}

void logics(){
	
	int a = 0;
	int b = 1;
	int c = 2;
	int d = 4;
	
	if(a) cout << "a" << endl;
	if(b) cout << "b" << endl;
	if(c) cout << "c" << endl;
	if(d) cout << "d" << endl;
	
	if(a>0 && b>0) cout << "ab" << endl;
	if(a>0 || d>0) cout << "ad" << endl;
	
}


void test(){
	double start_time = time(NULL); // measure running time
	
	
	//arrays();
	//byreference();
	//logics();
	//vectors();
	
	
	
	double end_time = 1.*( time(NULL));
	cout << "... runtime ~" << (end_time - start_time)/1. << " sec." << endl;
} // void test
