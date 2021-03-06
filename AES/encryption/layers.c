#include "layers.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>


static const unsigned char SBox[16][16] =  {
    {0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76},
    {0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0},
    {0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15},
    {0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75},
    {0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84},
    {0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf},
    {0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8},
    {0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2},
    {0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73},
    {0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb},
    {0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79},
    {0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08},
    {0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a},
    {0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e},
    {0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf},
    {0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16}
};

static const unsigned char MixColumnMat[4][4] =  {
	{2,3,1,1},
	{1,2,3,1},
	{1,1,2,3},
	{3,1,1,2}
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
				unsigned char judge = tmp / 128;
				tmp *= 2;
				if (judge > 0){
					tmp ^= num;
				}
			}
			ans ^= tmp;
		}
		b /= 2;
		i++;
	}
	return ans;
}

void generateRoundKey(unsigned char RoundKey[11][4][4], const unsigned char Key[4][4]){
	memcpy(RoundKey[0],Key,sizeof(unsigned char)*16);
	unsigned char r = 1;
	for (int i=1; i<11; i++){
		if(i!=1) 
			r *= 2;
		if(i>8)
			r ^= 0x1b;
		
		unsigned char a = SBox[RoundKey[i-1][1][3]/16][RoundKey[i-1][1][3]%16];
		unsigned char b = SBox[RoundKey[i-1][2][3]/16][RoundKey[i-1][2][3]%16];
		unsigned char c = SBox[RoundKey[i-1][3][3]/16][RoundKey[i-1][3][3]%16];
		unsigned char d = SBox[RoundKey[i-1][0][3]/16][RoundKey[i-1][0][3]%16];
		a = a^r;
		RoundKey[i][0][0] = RoundKey[i-1][0][0]^a;
		RoundKey[i][1][0] = RoundKey[i-1][1][0]^b;
		RoundKey[i][2][0] = RoundKey[i-1][2][0]^c;
		RoundKey[i][3][0] = RoundKey[i-1][3][0]^d;
		for (int j=0; j<4; j++){
			for (int k=1; k<4; k++){
				RoundKey[i][j][k] = RoundKey[i-1][j][k]^RoundKey[i][j][k-1];
			}
		}
	}
	/*
	for (int i=0; i<4; i++){
		for(int j=0; j<4; j++)
		printf("%x\n",RoundKey[1][j][i]);
	}*/
}

void SubBytes(unsigned char a[4][4]){
	for (int i=0 ; i<4; i++ ){
		for (int j=0; j<4; j++){
			unsigned char m = a[i][j]/16;
			unsigned char n = a[i][j]%16;
			a[i][j] = SBox[m][n];
		}
		
	}
}

void ShiftRows(unsigned char a[4][4]){
	for (int i=1; i<4; i++){
		unsigned char tmp[4];
		memcpy(tmp,a[i],sizeof(unsigned char)*4);
		for (int j=0; j<4; j++){
			a[i][j] = tmp[(j+i)%4];
		}
	}
}

void MixColumns(unsigned char a[4][4]){
	unsigned char ans[4][4];
	unsigned char value;
	for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            value = 0;
            for (int k = 0; k < 4; k++)
            {
                value ^= multiply(MixColumnMat[i][k],a[k][j]);
            }
            ans[i][j] = value;
        }
    }
    memcpy(a,ans,sizeof(unsigned char)*16);
}

void AddRoundKey(unsigned char a[4][4], unsigned char RoundKey[4][4]){
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			a[i][j] ^= RoundKey[i][j];
		}
	}
}

/*
int main(){
	unsigned char a[4][4] = {
		{0xe0, 0xd4,0xd4,0xd4},
		{0xb4, 0xbf,0xbf,0xbf},
		{0x52, 0x5d,0x5d,0x5d},
		{0xae, 0x30,0x30,0x30}};


	MixColumns(a);
	for (int i=0; i<4; i++){
		for (int j=0; j<4; j++){
			printf("%x\n",a[i][j]);
		}
	}
	
}*/
