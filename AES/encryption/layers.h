#ifndef LAYERS_H
#define LAYERS_H


void generateRoundKey(unsigned char RoundKey[11][4][4], const unsigned char Key[4][4]);

void SubBytes(unsigned char a[4][4]);

void ShiftRows(unsigned char a[4][4]);

void MixColumns(unsigned char a[4][4]);

void AddRoundKey(unsigned char a[4][4], unsigned char RoundKey[4][4]);

#endif
