#include <stdio.h> 
#include <stdlib.h> 
#include <iostream>
#include <gmp.h>
using namespace std;

void update(mpz_t m, mpz_t n, mpz_t q){
//Used to calculated "<m,n> <- <n,m-q*n>".
	mpz_t t;
	mpz_init(t);
	mpz_set(t,n);
	mpz_mul(n,q,t);
	mpz_sub(n,m,n);
	mpz_set(m,t);
	mpz_clear(t);
}


void EEA(mpz_t gcd, mpz_t a, mpz_t b){
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
	cout<<"s1:"<<s1<<endl;
	mpz_clear(s1);
	mpz_clear(t0);
	mpz_clear(t1);
	mpz_clear(q);
}


void rand_num(mpz_t num1,mpz_t num2){
//Used to genrate 2 random 4096-bit numbers.
	unsigned long seed;
	seed = time(NULL);
	gmp_randstate_t rstate;
	gmp_randinit_default(rstate);
	gmp_randseed_ui(rstate,seed);
	mpz_urandomb(num1,rstate,4096);
	mpz_urandomb(num2,rstate,4096);
}

int main () { 
	mpz_t num1, num2;
	mpz_init(num1);
	mpz_init(num2);
	mpz_set_ui(num1, 5993);
	mpz_set_ui(num2,31846);
	mpz_t result1,result2;
	mpz_init(result1);
	mpz_init(result2);
	mpz_gcd(result1,num1,num2);
	cout<<"The first generated random number is "<<num1<<"."<<endl;
	cout<<"The second generated random number is "<<num2<<"."<<endl;
	cout<<"GMP function result for gcd is "<<result1<<"."<<endl;
	EEA(result2,num1,num2);
	cout<<"Implemented extended Euclidean algorithm result for gcd is "<<result2<<"."<<endl;
	if(mpz_cmp(result1,result2) == 0)
		cout<<"The results are the same."<<endl;
	else
		cout<<"The results are different."<<endl;
	mpz_clear(num1);
	mpz_clear(num2);
	mpz_clear(result1);
	mpz_clear(result2);
} 



