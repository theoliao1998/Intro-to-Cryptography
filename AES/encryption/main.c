#include <stdio.h>
#include <string.h>
#include "layers.h"

static unsigned char plaintext[4][4] = 
{
	0,1,2,3,
	4,5,6,7,
	8,9,0,1,
	2,3,4,5
};

static unsigned char key[4][4] = 
{
   0,1,2,3,
   4,5,6,7,
   8,9,0,1,
   2,3,4,5
};

void AES(unsigned char ciphertext[4][4], unsigned char plaintext[4][4], unsigned char key[4][4]){
	unsigned char RoundKey[11][4][4];
	unsigned char a[4][4];
	memcpy(a,plaintext,sizeof(unsigned char)*16);
	generateRoundKey(RoundKey, key);
	AddRoundKey(a,RoundKey[0]);
	for(int i=1; i<=9; i++){
		SubBytes(a);
		ShiftRows(a);
		MixColumns(a);
		AddRoundKey(a,RoundKey[i]);
	}
	SubBytes(a);
	ShiftRows(a);
	MixColumns(a);
	AddRoundKey(a,RoundKey[10]);
	memcpy(ciphertext,a,sizeof(unsigned char)*16);
}


int main()
{
	unsigned char ciphertext[4][4];
	AES(ciphertext, plaintext, key);
	
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			printf("%x ",ciphertext[i][j]);
		}
		printf("\n");
	}
}

