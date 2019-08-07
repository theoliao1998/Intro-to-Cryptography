#include <math.h>
#include <stdio.h>

static const unsigned char M[8][8] = {
	1, 0, 0, 0, 1, 1, 1, 1,
	1, 1, 0, 0, 0, 1, 1, 1,
	1, 1, 1, 0, 0, 0, 1, 1,
	1, 1, 1, 1, 0, 0, 0, 1,
	1, 1, 1, 1, 1, 0, 0, 0,
	0, 1, 1, 1, 1, 1, 0, 0,
	0, 0, 1, 1, 1, 1, 1, 0,
	0, 0, 0, 1, 1, 1, 1, 1
};

static unsigned char multiply(unsigned char a, unsigned char b){
	unsigned char num = 0x1b;
	int i = 0;
	unsigned char ans = 0;
	unsigned char tmp;
	while (b > 0){
		if(b%2 != 0){
			tmp = a;
			for(int j=0; j<i; j++){
				unsigned char judge = tmp >> 7;
				tmp = tmp << 1;
				if (judge > 0){
					tmp ^= num;
				}
			}
			ans ^= tmp;
		}
		b = b >> 1;
		i++;
	}
	return ans;
}

unsigned int index_of_max(int a){
	int t = 1;
	int count = 0;
	for(unsigned int i =0; i<sizeof(int)*8; ++i){
		if(a&t) count = i;
		t = t << 1;
	}
	return count;
}

static int divide(int a, int b, int& r){
	int a_i = index_of_max(a);
	int b_i = index_of_max(b);
	if(a_i < b_i){
		r = a;
		return 0;
	}
	int c = a_i - b_i;
	a = a^(b<<c);
	return (1<<c) | divide(a,b,r);
	
}

static void update(int & a, int & b, int q){
	int v = 0;
	int tmp = b;
	int t = 1;
	for(unsigned int i=0; i<sizeof(int)*8;++i){
		if(q & t)
			v ^= ((b<<i));
		t = t << 1;
	}
	b = a^v;
	a = tmp;
}


static int EEA_inverse(unsigned char a){
//The implementation of extended Euclidean Algorithm.
	int r0,r1,s0,s1,t0,t1,q,r;
	r0 = 0x11b;
	r1 = a;
	s0 = 0;
	s1 = 1;
	t0 = 1;
	t1 = 0;
	while(r1!=1 && r1!=0){
		r = 0;
		q = divide(r0, r1, r);
		r0 = r1;
		r1 = r;
		update(s0, s1, q);
		update(t0, t1, q);
		//printf("%d\n",s1);
	}
	if (r1 == 0) return 0;
	return s1;
}

void generateSBox(unsigned char S[16][16]){
	unsigned char a,b;
	unsigned char b_mat[8];
	for(int i=0; i<16; i++){
		for(int j=0; j<16; j++){
			a = 16*i + j;
			
			//printf("%x\n",a);
			b = EEA_inverse(a);
			for(int k=0; k<8; k++){
				b_mat[k] = b%2;
				b /= 2;
			}

			unsigned char c = 0;
			unsigned char d = 1;
			for (int m = 0; m < 8; m++){
				unsigned char value = 0;
				for (int n = 0; n < 8; n++){
					value ^= multiply(M[m][n],b_mat[n]);
				}
				c ^= (value%2)*d;
				d *= 2;
			}
			c ^= 99;

			S[i][j] = c;
		}
	}
}

int main(){
	unsigned char S[16][16];
	//printf("%d",EEA_inverse(0));
	
	generateSBox(S);
	printf("S-Box\n");
	for (int i=0; i<16; i++){
		for (int j=0; j<16; j++){
			printf("%d ",S[i][j]);
		}
		printf("\n");
	}
}
