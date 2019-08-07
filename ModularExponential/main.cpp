#include <stdio.h> 
#include <stdlib.h> 
#include <iostream>
#include <time.h>
#include <gmp.h>
using namespace std;


void ModularExponential(mpz_t power,const mpz_t m,const mpz_t d,const mpz_t n){
	mpz_set_ui(power,1);	
	auto k = mpz_sizeinbase(d, 2);

	for (int i = k-1; i>=0; i--){
		mpz_mul(power,power,power);
		mpz_mod(power,power,n);
		if (mpz_tstbit(d,i)){
			mpz_mul(power,m,power);
			mpz_mod(power,power,n);
		}
	}
	cout<<endl;
}

void RewritedModularExponential(mpz_t power,const mpz_t alpha,const mpz_t a, const mpz_t beta, const mpz_t b, const mpz_t n){
	mpz_t c,prod,d,e,power1,power2;
	mpz_init(c);
	mpz_init(d);
	mpz_init(e);
	mpz_init(prod);
	mpz_init(power1);
	mpz_init(power2);
	mpz_set_ui(power1,1);
	mpz_set_ui(power2,1);
	mpz_mul(prod,alpha,beta);
	if(mpz_cmp(a,b) >= 0){
		mpz_set(c,b);
		mpz_sub(d,a,b);
		mpz_set(e,alpha);
	}
	else{
		mpz_set(c,a);
		mpz_sub(d,b,a);
		mpz_set(e,beta);
	}
	auto k1 = mpz_sizeinbase(c, 2);
	auto k2 = mpz_sizeinbase(d, 2);
	for(int i = k1-1; i>=0; i--){
		mpz_mul(power1,power1,power1);
		mpz_mod(power1,power1,n);
		if (mpz_tstbit(c,i) == 1){
			mpz_mul(power1,prod,power1);
			mpz_mod(power1,power1,n);
		} 
	}
	for(int i = k2-1; i>=0; i--){
		mpz_mul(power2,power2,power2);
		mpz_mod(power2,power2,n);
		if (mpz_tstbit(d,i) == 1){
			mpz_mul(power2,e,power2);
			mpz_mod(power2,power2,n);
		} 
	}
	mpz_mul(power,power1,power2);
	mpz_mod(power,power,n);
	mpz_clears(c,prod,d,e,power1,power2,NULL);
}

const int L = 50;

void rand_num(mpz_t num1,mpz_t num2,mpz_t num3,mpz_t num4,mpz_t num5){
//Used to genrate 5 random L-bit numbers.
	unsigned long seed;
	seed = time(NULL);
	gmp_randstate_t rstate;
	gmp_randinit_default(rstate);
	gmp_randseed_ui(rstate,seed);
	mpz_urandomb(num1,rstate,L);
	mpz_urandomb(num2,rstate,L);
	mpz_urandomb(num3,rstate,L);
	mpz_urandomb(num4,rstate,L);
	do{
		mpz_urandomb(num5,rstate,L);}
	while (mpz_cmp_ui(num5,0) == 0);
		
}


int main(){
	mpz_t alpha, a, beta, b, n, power1, power2, power3;
	mpz_inits(alpha, a, beta, b, n, power1, power2, power3, 0);
	rand_num(alpha, a, beta, b, n);
	clock_t start = clock();
	ModularExponential(power1, alpha, a, n);
	ModularExponential(power2, beta, b, n);
	mpz_mul(power3,power1,power2);
	mpz_mod(power3,power3,n);
	clock_t end = clock();
	cout<<"To calculate "<<alpha<<"^"<<a<<"*"<<beta<<"^"<<b<<" mod "<<n<<" using the original strategy, it costs "<<end-start<<"ms. The result is "<<power3<<"."<<endl;
	start = clock();
	RewritedModularExponential(power3, alpha, a, beta, b, n);
	end = clock();
	cout<<"To calculate "<<alpha<<"^"<<a<<"*"<<beta<<"^"<<b<<" mod "<<n<<" using the rewrited strategy, it costs "<<end-start<<"ms. The result is "<<power3<<"."<<endl;

}
