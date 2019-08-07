#include <stdio.h> 
#include <stdlib.h> 
#include <iostream>
#include <gmp.h>
using namespace std;

struct key_p_q{
	mpz_t p;
	mpz_t q;
	mpz_t n;
};

static void rand_prime(mpz_t num1,mpz_t num2, int level){
	int length = level == 80 ? 1024:
				 level == 112 ? 2048:
				 level == 128 ? 3072:
				 level == 192 ? 7680:
				 level == 256 ? 15360 : 8;
	unsigned long seed;
	seed = time(NULL);
	gmp_randstate_t rstate;
	gmp_randinit_default(rstate);
	gmp_randseed_ui(rstate,seed);
	mpz_rrandomb(num1,rstate,length/2);
	mpz_nextprime(num1,num1);
	mpz_rrandomb(num2,rstate,length/2);
	mpz_nextprime(num2,num2);	
}

key_p_q generate(int level){
	key_p_q key;
	mpz_init(key.p);
	mpz_init(key.q);
	rand_prime(key.p,key.q,level);
	mpz_init(key.n);
	mpz_mul(key.n,key.p,key.q);
	return key;
}

void encrypt(mpz_t c, mpz_t e, mpz_t m, mpz_t n){
	mpz_powm(c, m, e, n);
}

void decrypt(mpz_t m, mpz_t d, mpz_t c, mpz_t n){
	mpz_powm(m, c, d, n);
}

void RSA(int level, mpz_t m){
	cout<<"The security level is "<<level<<"."<<endl;
	cout<<"The message is "<<m<<"."<<endl;
	
	//generate p,q, and n
	key_p_q k_p_q = generate(level);

	//calculate \varphi(n)
	mpz_t phi,p_1,q_1;
	mpz_init(phi);
	mpz_init(p_1);
	mpz_init(q_1);
	mpz_sub_ui(p_1,k_p_q.p,1);
	mpz_sub_ui(q_1,k_p_q.q,1);
	mpz_mul(phi,p_1,q_1);
	

	//generate e and d
	mpz_t e,d,gcd;
	mpz_init(e);
	mpz_init(d);
	mpz_init(gcd);	
	unsigned long seed;
	seed = time(NULL);
	gmp_randstate_t rstate;
	gmp_randinit_default(rstate);
	gmp_randseed_ui(rstate,seed);
	while(1){
		mpz_urandomm(e,rstate,phi);
		mpz_gcd(gcd,e,phi);
		if(mpz_cmp_ui(e,0) > 0 && mpz_cmp_ui(gcd,1) == 0)
			break;
	}
	mpz_invert(d, e, phi);
	
	cout<<"The public keys are: n = "<<k_p_q.n<<", e = "<<e<<"."<<endl;
	cout<<"The private key is: d = "<<d<<"."<<endl;

	//encrypt
	mpz_t c;
	mpz_init(c);
	encrypt(c,e,m,k_p_q.n);
	
	cout<<"The ciphertext is "<<c<<"."<<endl;
	
	//decrypt
	mpz_t M;
	mpz_init(M);
	decrypt(M,d,c,k_p_q.n);
	
	cout<<"After decryption, the text becomes "<<M<<"."<<endl;
	
	mpz_clears(k_p_q.p,k_p_q.q,k_p_q.n,phi,p_1,q_1,e,d,gcd,c,M,NULL);
}

int main(){
	mpz_t m;
	mpz_init(m);
	mpz_set_ui(m,12343);
	RSA(80,m);
}
