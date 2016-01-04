
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



void test(){
	
	
	//arrays();
	byreference();
	
	
	cout << "done" << endl;
	
	
	
	
} // void test
