#include <stdio.h> 
#include <stdlib.h> 
#include <iostream>
#include <gmp.h>
using namespace std;

static void f(mpz_t ans, const mpz_t x, const mpz_t n){
	mpz_mul(ans,x,x);
	mpz_add_ui(ans,ans,1);
	mpz_mod(ans,ans,n);
}

static void update(mpz_t m, mpz_t n, mpz_t q){
//Used to calculated "<m,n> <- <n,m-q*n>".
	mpz_t t;
	mpz_init(t);
	mpz_set(t,n);
	mpz_mul(n,q,t);
	mpz_sub(n,m,n);
	mpz_set(m,t);
	mpz_clear(t);
}

static void gcd_EEA(mpz_t gcd, const mpz_t a, const mpz_t b){
//The implementation of extended Euclidean Algorithm.
	mpz_t r0,r1,s0,s1,t0,t1,q;
	mpz_init(r0);
	mpz_init(r1);
	mpz_init(s0);
	mpz_init_set_str(s1, "1", 10);
	mpz_init_set_str(t0, "1", 10);
	mpz_init(t1);
	mpz_init(q);
	mpz_set(r0,b);
	mpz_set(r1,a);

	while(mpz_sgn(r0)!=0){
		mpz_tdiv_q(q,r1,r0);
		update(r1, r0, q);
		update(s1, s0, q);
		update(t1, t0, q);
	}
	mpz_set(gcd,r1);
	mpz_clear(r0);
	mpz_clear(r1);
	mpz_clear(s0);
	mpz_clear(s1);
	mpz_clear(t0);
	mpz_clear(t1);
	mpz_clear(q);
}

bool pollard_rho(mpz_t d, const mpz_t n){
	mpz_t a,b,a_b;
	mpz_init(a);
	mpz_init(b);
	mpz_init(a_b);
	mpz_set_ui(a,2);
	mpz_set_ui(b,2);
	mpz_set_ui(d,1);
	while(mpz_cmp_ui(d,1) == 0){
		f(a,a,n);
		f(b,b,n);
		f(b,b,n);
		if(mpz_cmp(a,b) >= 0)
			mpz_sub(a_b,a,b);
		else
			mpz_sub(a_b,b,a);
		gcd_EEA(d, a_b, n);
	}
	if(mpz_cmp(d,n) == 0)
		return 0;
	else 
		return 1;
}


int main(){
	unsigned long int number;
	cout<<"Please input an unsigned number."<<endl;
	cin>>number; 
	mpz_t n,d;
	mpz_init(n);
	mpz_init(d);
	mpz_set_ui(n,number);
	if (pollard_rho(d,n))
		cout<<d<<" is a nontrivial factor of "<<n<<"."<<endl;
	else
		cout<<"Failure."<<endl;
}
