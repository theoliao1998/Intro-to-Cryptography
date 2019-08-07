#include <stdio.h> 
#include <string.h>
#include <stdlib.h>
#include <math.h> 
#include <limits>
#include <iostream>
#include <gmp.h>
#include <vector>
#include <map>
#include <time.h>

typedef unsigned long int u_type;

using namespace std;

u_type nextPrime(u_type p){
	while(1){
		p += 1;
		u_type bound = sqrt(p);
		bool prime = 1;
		for(u_type i = 2; i<=bound; i++){
			if(p%i == 0){
				prime = 0;
				break;
			}
		}
		if(prime==1){
			break;
		}
	}
	return p;
}

static void getParameters(const mpz_t N, u_type & M, u_type & F, double & T){
	//parameters selected according to Silverman's paper
	auto length = mpz_sizeinbase(N,10);
	if(length <= 24){
		F = 100;
		M = 5000;
		T = 1.5;
	}
	else if(length <= 30){
		F = 200;
		M = 25000;
		T = 1.5;
	}
	else if(length <= 36){
		F = 400;
		M = 25000;
		T = 1.75;
	}
	else if(length <= 42){
		F = 900;
		M = 50000;
		T = 2.0;
	}
	else if(length <= 48){
		F = 1200;
		M = 100000;
		T = 2.0;
	}
	else if(length <= 54){
		F = 2000;
		M = 250000;
		T = 2.2;
	}
	else if(length <= 60){
		F = 3000;
		M = 350000;
		T = 2.4;
	}
	else if(length <= 66){
		F = 4500;
		M = 500000;
		T = 2.6;
	}
	
}

static void getFB(mpz_t N, const u_type F,  u_type FB[]){
	u_type p = 2;
	u_type i = 0;
	while(i<F){
		if(mpz_kronecker_ui (N, p) == 1){
			FB[i] = p;
			i++;
		}
		p = nextPrime(p);
	}
	
}


u_type getMultiplier(mpz_t N, u_type F){
	mpz_t kN;
	u_type p;
	u_type res=1;
	double f,g,fmax;
	fmax = (numeric_limits<double>::min)();
	mpz_init(kN);
	for (u_type k=1; k<100; k++){
		mpz_mul_ui(kN,N,k);
		f = 0;
		u_type* FB;
		FB = new u_type[F];
		getFB(kN,F,FB);
		for(u_type i = 0; i<F; i++){
			p = FB[i];
			if(p==2)
				g = mpz_fdiv_ui (kN, 8) == 1 ? 2.0 : 0;
			else
				g = k % p == 0 ? (1.0/p) : (2.0/p);
			f += g * log((double)p);
		}
		f -= 0.5 * log((double)k);
		if(f>fmax){
			fmax= f;
			res = k;
		}
	}
	mpz_clear (kN);
	return res;
}

static long int getTestValue(u_type M,mpz_t kN, u_type pmax,u_type T){
	mpf_t res,KN,pmax_T;
	mpf_inits(res,KN,pmax_T,0);
	mpf_set_z(KN,kN);
	mpf_set_d(pmax_T,pow(pmax,T));
	mpf_div_ui(res,KN,2);
	mpf_sqrt(res,res);
	mpf_mul_ui(res,res,M);
	mpf_div(res,res,pmax_T);
	double ans = mpf_get_d(res);
	ans = log(ans);
	mpf_clears(res,KN,pmax_T,NULL);
	return (long int) ans;
}

static bool goodD(const mpz_t D, const mpz_t kN){
	return mpz_fdiv_ui(D,4) == 3 && mpz_kronecker(D,kN) == 1;
}

static void selectCoefficients(mpz_t A, mpz_t B, mpz_t C, mpz_t D, const mpz_t D0, const mpz_t kN, bool firstround){

	mpz_t t;
	mpz_init(t);
	u_type distance = 0;
	if(!firstround){
		mpz_sub(t,D,D0);
		mpz_abs(t,t);
		distance = mpz_get_ui(t);
	}
	
	for(u_type i = distance; 1 ; i++){
		mpz_add_ui(t,D0,i);
		if(goodD(t,kN) && (firstround || mpz_cmp(t,D)!=0)){
			mpz_add_ui(D,D0,i);
			break;
		}
		mpz_sub_ui(t,D0,i);
		if( mpz_cmp_ui(D0,i)>=0 && goodD(t,kN) && (firstround || mpz_cmp(t,D)!=0)){
			mpz_sub_ui(D,D0,i);
			break;
		}
	}
	
	mpz_pow_ui(A,D,2);
	
	mpz_t h0,h1,h2;
	mpz_inits(h0,h1,h2,0);
	mpz_sub_ui(h0,D,3);
	mpz_fdiv_q_ui(h0,h0,4);
	mpz_powm(h0,kN,h0,D);
	mpz_add_ui(h1,D,1);
	mpz_fdiv_q_ui(h1,h1,4);
	mpz_powm(h1,kN,h1,D);
	mpz_powm_ui(t,h1,2,D);
	mpz_sub(t,kN,t);
	mpz_fdiv_q(t,t,D);
	mpz_mul_ui(h2,h1,2);
	mpz_invert(h2,h2,D);
	mpz_mul(h2,h2,t);
	mpz_mod(h2,h2,D);
	mpz_mul(B,h2,D);
	mpz_add(B,h1,B);
	mpz_mod(B,B,A);
	if(mpz_odd_p(B) == 0)
		mpz_sub(B,A,B);
	mpz_pow_ui(t,B,2);
	mpz_sub(t,t,kN);
	mpz_mul_ui(C,A,4);
	mpz_fdiv_q(C,t,C);
	
	mpz_clears(h0,h1,h2,t,NULL);
}

static u_type mpz_sqrtm(mpz_t x, const mpz_t n, const mpz_t p){
	//implement Tonelli-Shanks algorithm, p is prime
	if(mpz_cmp_ui(p,2)==0){
		if(mpz_odd_p(x) == 0){
			mpz_set_ui(x,0);
		}
		else{
			mpz_set_ui(x,1);
		}
		return 1;
	}
	if(mpz_legendre(n, p) != 1)
		//should not happen for FB selected
		return 0;
	
	if(mpz_fdiv_ui(p,4)==3){
		mpz_add_ui(x,p,1);
		mpz_fdiv_q_ui(x,x,4);
		mpz_powm(x,n,x,p);
		return 1;
	}
	mpz_t q,s,z,c,r,t,m,b,tmp;
	mpz_inits(q,s,z,c,r,t,m,b,tmp,0);
	mpz_sub_ui(q,p,1);
	while(mpz_odd_p(q) == 0){
		mpz_add_ui(s,s,1);
		mpz_fdiv_q_ui(q,q,2);
	}
	mpz_set_ui(z,2);
	while(mpz_legendre(z, p) != -1) 
		mpz_add_ui(z, z, 1);            
	mpz_powm(c,z,q,p);
	mpz_add_ui(r,q,1);
	mpz_fdiv_q_ui(r,r,2);
	mpz_powm(r,n,r,p);
	mpz_powm(t,n,q,p);
	mpz_set(m,s);
	while(1){
		if(mpz_cmp_ui(t,1) == 0){
			mpz_set(x,r);
			break;
		}
		else{
			u_type i;
			mpz_set(tmp,t);
			for(i=1;mpz_cmp_ui(m,i)>0;i++){
				mpz_powm_ui(tmp,tmp,2,p);
				if(mpz_cmp_ui(tmp,1)==0)
					break;
			}
			mpz_sub_ui(b,m,(i+1));
			u_type v = mpz_get_ui(b);
			mpz_ui_pow_ui(b,2,v);
			mpz_powm(b,c,b,p);
			mpz_mul(r,r,b);
			mpz_mod(r,r,p);
			mpz_powm_ui(c,b,2,p);
			mpz_mul(t,t,c);
			mpz_mod(t,t,p);
			mpz_set_ui(m,i);
		}
	}
	mpz_clears(q,s,z,c,r,t,m,b,tmp,NULL);
	return 1;
}

static void sieve(const mpz_t A,const mpz_t B, const mpz_t kN, const u_type FB[], u_type F, long int SieveArray[],u_type M){
	mpz_t p,t,r1,r2;
	mpz_inits(p,t,r1,r2,0);
	u_type R1,R2;
	for(u_type i=0;i<F;i++){
		mpz_set_ui(p,FB[i]);
		mpz_mul_ui(t,A,2);
		mpz_invert(t,t,p);
		mpz_sqrtm(r1,kN,p);
		mpz_sub(r2,p,r1);
		mpz_sub(r1,r1,B);
		mpz_sub(r2,r2,B);
		mpz_mul(r1,r1,t);
		mpz_mul(r2,r2,t);
		mpz_mod(r1,r1,p);
		mpz_mod(r2,r2,p);
		R1 = mpz_get_ui(r1);
		R2 = mpz_get_ui(r2);

		u_type P = mpz_get_ui(p);
		long int v = log(1.0*P);
		if(R1 <= M)
			SieveArray[M+R1] += v;
		if(R2 <= M)
			SieveArray[M+R2] += v;
		long int j = M+R1+P;
		
		while(j<(int)(2*M+1)){
			SieveArray[j] += v;
			j += P;
		}
		j = M+R1-P;
		while(j>=0){
			SieveArray[j] += v;
			j -= P;
		}
		
		j = M+R2+P;
		while(j<(int)(2*M+1)){
			SieveArray[j] += v;
			j += P;
		}
		j = M+R2-P;
		while(j>=0){
			SieveArray[j] += v;
			j -= P;
		}
		
	}
	
	mpz_clears(p,t,r1,r2,NULL);
}

static void swapRows(u_type** Matrix,u_type** operations, u_type F,u_type size, u_type i, u_type j){
	if(i==j)
		return;
	else{
		u_type tmp1[F];
		u_type tmp2[size];
		
		memcpy(tmp1,Matrix[i],sizeof(tmp1));
		memcpy(tmp2,operations[i],sizeof(tmp2));
		memcpy(Matrix[i],Matrix[j],sizeof(tmp1));
		memcpy(operations[i],operations[j],sizeof(tmp2));
		memcpy(Matrix[j],tmp1,sizeof(tmp1));
		memcpy(operations[j],tmp2,sizeof(tmp2));
	}
}

static void addRows(u_type** Matrix,u_type** operations, u_type F,u_type size,u_type dest, u_type a, u_type b){
	for(u_type i=0; i<F; i++){
		Matrix[dest][i] = Matrix[a][i]+Matrix[b][i];
	}
	for(u_type i=0; i<size; i++){
		operations[dest][i] = operations[a][i]+operations[b][i];
	}
	
}

static long int checkEven(u_type** Matrix,u_type** operations, u_type F,u_type size,u_type i, bool LisOne){
	long int res = -1;
	for(u_type j=i; j<size; j++){
		u_type k;
		for(k=0; k<F && (Matrix[j][k]%2)==0;k++);
		if(k==F){
			if(LisOne) return j;
			u_type op_num = 0;
			for(u_type l=0; l<size; l++){
				op_num += operations[j][l];
			}
			if((op_num%2)==0){
				res = j;				
				break;
			}
		}
	}
	return res;
}


static long int GaussianEliminate(u_type** exp_vectors,u_type** operations, u_type F, u_type size, bool LisOne){
	long int res = -1;
	
	for(u_type i=0; i<size; i++){
		u_type j;
		for(j=i; j<size; j++){
			if((exp_vectors[j][i]%2) == 1){
				swapRows(exp_vectors,operations,F,size, j, i);
				break;
			}
		}
		
		if(j<size){
			for(j=i+1; j<size; j++){
				if((exp_vectors[j][i]%2) == 1){
					addRows(exp_vectors,operations,F,size,j, j, i);
				}
			}
		}
		res = checkEven(exp_vectors,operations,F,size,i,LisOne);
		
		if(res>=0){
			return res;
		}
	}
	return res;
}


int main() { 
	clock_t start = clock();
	clock_t end;
	clock_t TimeLimit = 15 * CLOCKS_PER_SEC; //time limit, should be modified according to the length of N
	
	mpz_t N,kN;
	u_type M,F,k;
	M = 0;
	F = 0;
	k = 0;
	double T;
	T = 0;
	//initialize N
	mpz_init(N);
	unsigned long seed;
	seed = time(NULL);
	gmp_randstate_t rstate;
	gmp_randinit_default(rstate);
	gmp_randseed_ui(rstate,seed);
	mpz_urandomb(N,rstate,24);
	
	
	//eliminate small prime factors
	for(u_type i=2; i<=5;i++){
		while(mpz_divisible_ui_p(N,i)!=0){
			mpz_divexact_ui(N,N,i);
		}
	}

	
	cout<<"N = "<<N<<endl;
	
	
	
	//Generate M,F,T
	getParameters(N,M,F,T);
	
	//Determine multiplier k
	k = getMultiplier(N, F);
	
	mpz_init(kN);
	mpz_mul_ui(kN,N,k);
	
	//Generate FB
	u_type * FB;
	FB = new u_type[F];
	getFB(kN,F,FB);
    
    //Compute test value
    long int testvalue = getTestValue(M,kN,FB[F-1],T);
    
    
    //initialize coeffients A,B,C where Q(x)=Ax^2+Bx+C
    
    mpz_t A,B,C,D;
    mpz_inits(A,B,C,D,0);
    
    mpf_t d;
	mpf_init(d);
	mpf_set_z(d,kN);
	mpf_div_ui(d,d,2);
	mpf_sqrt(d,d);
	mpf_div_ui(d,d,M);
	mpf_sqrt(d,d);
	
	mpz_t D0;
	mpz_init(D0);
	mpz_set_f(D0,d);
	mpf_clear(d);
    bool firstround = 1;
    
	mpz_t* H;
	H = new mpz_t[2*M+1];//save H for each i
	for(u_type i=0; i< 2*M+1; i++){
		mpz_init(H[i]);
	}
	
	 
	u_type** V;
	
	mpz_t* residues;
	residues = new mpz_t[2*M+1];
	
	V = new u_type*[2*M+1];// save base vector
	//the first F values correspond to exponent vectors, and the last is 
    
    for(u_type i=0; i< 2*M+1; i++){
		V[i] = new u_type[F];
		mpz_init(residues[i]);
		for(u_type j=0; j< F; j++){
			V[i][j] = 0;
		}
	}
    
    bool EnoughFactorized = 0;
    u_type FT = 0; //number of factoraizations
    u_type FF = 0; //number of residues fully factored 
   
    
    bool* factorized;
    factorized = new bool[2*M+1];
    for(u_type i=0; i<2*M+1; i++){
		factorized[i] = 0;
	}
    
    mpz_t H_tmp,P1,P2;
	mpz_inits(H_tmp,P1,P2,0);
	mpz_t Q_x;
	mpz_init(Q_x);
	mpz_t tmp;
	mpz_init(tmp);
	mpz_t p; //prime p_i in FB
	mpz_init(p);
    
    while(!EnoughFactorized){
		end = clock();
		if((end-start)>TimeLimit) break;

		//Select coefficients
		selectCoefficients(A,B,C,D,D0,kN,firstround);
		firstround = 0;		
		//initialize sieve array
		
		long int SieveArray[2*M+1] = {0};
		sieve(A,B,kN,FB,F,SieveArray,M);
		
		for(u_type i = 0; i < (2*M+1); i++){
			if(SieveArray[i]>testvalue){
				//calculate H and Q_x
				mpz_mul_ui(Q_x,A,i-M);
				mpz_mul_ui(Q_x,Q_x,2);
				mpz_add(Q_x,Q_x,B);
				mpz_mul_ui(tmp,D,2);
				mpz_invert(tmp,tmp,kN);
				mpz_mul(Q_x,Q_x,tmp);
				mpz_mod(Q_x,Q_x,kN);
				mpz_set(H[i],Q_x);
				mpz_powm_ui(Q_x, Q_x, 2,kN);
				
				mpz_set(H_tmp,Q_x);

				//factorize through division
				for (u_type j = 0; j<F;j++){
					mpz_set_ui(p,FB[j]);
					mpz_mod(tmp, Q_x, p);
					bool isBase = (mpz_cmp_ui(tmp, 0) == 0);
					while(isBase)
					{
						mpz_fdiv_q(Q_x,Q_x,p);
						V[i][j]++;
						mpz_mod(tmp, Q_x, p);
						isBase = (mpz_cmp_ui(tmp, 0) == 0);
					}
				}
				
				mpz_set(residues[i],Q_x);
				
						
				
				if (mpz_cmp_ui(residues[i], 1) == 0) // fully factorized
				{	
					factorized[i] = 1;
					FF++;
					FT++;
				}
				else // not fully factorized
				{
					if(mpz_cmp(residues[i],H_tmp)!=0){
						factorized[i] = 1;
						FT++;
					}
				}
				
			}
		}
		
		if(((1.0*FT)/(F+FT-FF)) >= 0.96){
		
		map<string, vector<u_type> > mapL;
		map<string, vector<u_type>>::iterator it;
		string L;
	
		u_type i;
		for(i = 0; i<2*M+1; i++){
			if(!factorized[i]) continue;
		L = mpz_get_str(NULL,10,residues[i]);

		
	    it=mapL.find(L); //map L to a list of i s.
	    if(it==mapL.end()){
			vector<u_type> v;
			v.push_back(i);
			mapL[L] = v;
		}
		else{
			mapL[L].push_back(i);
			if(mapL[L].size() >= 2){
				
				u_type ** exp_vectors;
				u_type size = mapL[L].size();
				exp_vectors = new u_type*[size];
				for(u_type f=0; f<size; f++){
					exp_vectors[f] = new u_type[F];
				}
				u_type ** operations;
				operations=new u_type*[size];
				for(u_type f=0;f<size;f++){
					operations[f]=new u_type[size];
				}
				
				
				
				for(u_type j=0; j<size; j++){
					for(u_type k=0; k<size; k++){
						operations[j][k] = 0;
					}
				}
				
				
				for(u_type j=0; j<size; j++){
					operations[j][j] = 1;
				}
				
				for(u_type j=0; j<size;j++){
					for(u_type k=0; k<F; k++)
						exp_vectors[j][k] = V[mapL[L][j]][k];
				}

				
				//if L is not 1, then we need even numbers multiplied together to obtain L^2
				bool LisOne = (mpz_cmp_ui(residues[i],1)==0);
				
				//Gaussian elimination
				long int res = GaussianEliminate(exp_vectors,operations,F,mapL[L].size(),LisOne);
				
				
				//obatin P1 and P2 according to the operations in Gaussain elimination
				if(res>=0){
					mpz_set_ui(H_tmp,1);
					for(u_type j=0; j<size; j++){
						mpz_powm_ui(tmp,H[mapL[L][j]],operations[res][j],kN);
						
						mpz_mul(H_tmp,H_tmp,tmp);
						mpz_mod(H_tmp,H_tmp,kN);
					}
					mpz_set(P1,H_tmp);
					mpz_set_ui(P2,1);
					for(u_type j=0; j<F; j++){
						mpz_set_ui(tmp,FB[j]);
						if(exp_vectors[res][j]>0){
							
							mpz_powm_ui(tmp,tmp,exp_vectors[res][j]/2,kN);
							mpz_mul(P2,P2,tmp);
							mpz_mod(P2,P2,kN);
						}
					}
					u_type op_num = 0;
					for(u_type l=0; l<size; l++){
						op_num += operations[res][l];
					}
					mpz_powm_ui(tmp,residues[i],op_num/2,kN);
					mpz_mul(P2,P2,tmp);
					mpz_mod(P2,P2,kN);
					
					
					
					mpz_add(tmp,P1,P2);
					mpz_sub(H_tmp,P1,P2);
					mpz_gcd(tmp,tmp,N);
					mpz_gcd(H_tmp,H_tmp,N);
					//gcd(P1+P2,N) or gcd(P1-P2,N) is likely to give a nontrivious factor of N
					
					if(mpz_cmp_ui(tmp,1) !=0 &&  mpz_cmp(tmp,N) !=0){
						cout<<tmp<<" is a nontrivious factor of "<<N<<endl;
						mpz_powm_ui(P1,P1,2,kN);
						mpz_powm_ui(P2,P2,2,kN);
						
						//we expect P1^2 \equiv P2^2 \bmod kN
						cout<<"P1^2 = "<<P1<<" P2^2 = "<<P2<<endl;
						for(u_type f=0; f<size; f++){
							delete []exp_vectors[f];
							delete []operations[f];
						}
						delete []exp_vectors;
						delete []operations;
						break;
					}
					else if(mpz_cmp_ui(H_tmp,1)!=0 &&  mpz_cmp(H_tmp,N) !=0){
						cout<<H_tmp<<" is a nontrivious factor of "<<N<<endl;
						mpz_powm_ui(P1,P1,2,kN);
						mpz_powm_ui(P2,P2,2,kN);
						cout<<"P1^2: "<<P1<<" P2^2: "<<P2<<endl;
						for(u_type f=0; f<size; f++){
							delete []exp_vectors[f];
							delete []operations[f];
						}
						delete []exp_vectors;
						delete []operations;
						break;
					}
					
					
				}
				
				for(u_type f=0; f<size; f++){
					delete []exp_vectors[f];
					delete []operations[f];
				}
				delete []exp_vectors;
				delete []operations;
			}
			
		}
	}
	 if(i<2*M+1){
		EnoughFactorized = 1;
	 }
	
	it = mapL.begin();
    while(it != mapL.end()) {
        it->second.clear();
        it++;
    }
	mapL.clear();
	}
		
	}
	 
    if (!EnoughFactorized)
		cout<<"Fail to find a non-trivious factor within the time limit. It might be prime."<<endl;
    
    //free
	mpz_clears(N,kN,A,B,C,D,D0,tmp,H_tmp,P1,P2,Q_x,p,NULL);
	for(u_type i=0; i<2*M+1; i++){
		mpz_clears(H[i],residues[i],NULL);
		delete []V[i];
	}
	delete []H;
	delete []residues;
	delete []V;
	delete []FB;
	delete []factorized;
}



