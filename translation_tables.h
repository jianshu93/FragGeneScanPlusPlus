#ifndef __TRANSLATION_TABLES_H
#define __TRANSLATION_TABLES_H

#include <stdbool.h>

int TRANSLATION_TABLES[20] = { 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26 };
#define NUM_TRANSLATION_TABLES 20;

char TRANSLATION_TABLE_1[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    '*', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_1_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    '*', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_2[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    '*', 'S', '*', 'S',
    'M', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    'W', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_2_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', '*',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', '*',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_3[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'M', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'T', 'T', 'T', 'T',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    'W', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_3_RC[65] = {
    'F', 'V', 'T', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'T', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'T', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'T', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_4[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    'W', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_4_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_5[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'S', 'S', 'S', 'S',
    'M', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    'W', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_5_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_6[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    'Q', 'Y', 'Q', 'Y',
    'S', 'S', 'S', 'S',
    '*', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_6_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    'Q', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    '*', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    'Q', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_9[65] = {
    'N', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'S', 'S', 'S', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    'W', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_9_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    'W', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'N',
    'X'
};

char TRANSLATION_TABLE_10[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    'C', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_10_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    'C', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_11[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    '*', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_11_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    '*', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_12[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'S', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    '*', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_12_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'S', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    '*', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_13[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'G', 'S', 'G', 'S',
    'M', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    'W', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_13_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'G',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'G',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_14[65] = {
    'N', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'S', 'S', 'S', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    'Y', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    'W', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_14_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    'W', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'E', 'Q', 'N',
    'X'
};

char TRANSLATION_TABLE_15[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', 'Q', 'Y',
    'S', 'S', 'S', 'S',
    '*', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_15_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    'Q', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    '*', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_16[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', 'L', 'Y',
    'S', 'S', 'S', 'S',
    '*', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_16_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    'L', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    '*', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_21[65] = {
    'N', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'S', 'S', 'S', 'S',
    'M', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    'W', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_21_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'N',
    'X'
};

char TRANSLATION_TABLE_22[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', 'L', 'Y',
    '*', 'S', 'S', 'S',
    '*', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_22_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    'L', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    '*', 'G', 'R', 'R',
    '*', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_23[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    '*', 'C', 'W', 'C',
    '*', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_23_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    '*', 'V', 'L', 'I',
    '*', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_24[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'S', 'S', 'K', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    'W', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_24_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'K',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    'W', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_25[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'L', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    'G', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_25_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    'G', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

char TRANSLATION_TABLE_26[65] = {
    'K', 'N', 'K', 'N',
    'T', 'T', 'T', 'T',
    'R', 'S', 'R', 'S',
    'I', 'I', 'M', 'I',
    'Q', 'H', 'Q', 'H',
    'P', 'P', 'P', 'P',
    'R', 'R', 'R', 'R',
    'L', 'L', 'A', 'L',
    'E', 'D', 'E', 'D',
    'A', 'A', 'A', 'A',
    'G', 'G', 'G', 'G',
    'V', 'V', 'V', 'V',
    '*', 'Y', '*', 'Y',
    'S', 'S', 'S', 'S',
    '*', 'C', 'W', 'C',
    'L', 'F', 'L', 'F',
    'X'
};

char TRANSLATION_TABLE_26_RC[65] = {
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'A', 'M',
    'W', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'F', 'V', 'L', 'I',
    'C', 'G', 'R', 'S',
    'S', 'A', 'P', 'T',
    'Y', 'D', 'H', 'N',
    'L', 'V', 'L', 'I',
    '*', 'G', 'R', 'R',
    'S', 'A', 'P', 'T',
    '*', 'E', 'Q', 'K',
    'X'
};

bool parse_translation_tables(int id, char *out_table[], char *out_table_rc[]) {
    switch (id) {
        case 1: *out_table = TRANSLATION_TABLE_1;
                 *out_table_rc = TRANSLATION_TABLE_1_RC;
                 return true;
        case 2: *out_table = TRANSLATION_TABLE_2;
                 *out_table_rc = TRANSLATION_TABLE_2_RC;
                 return true;
        case 3: *out_table = TRANSLATION_TABLE_3;
                 *out_table_rc = TRANSLATION_TABLE_3_RC;
                 return true;
        case 4: *out_table = TRANSLATION_TABLE_4;
                 *out_table_rc = TRANSLATION_TABLE_4_RC;
                 return true;
        case 5: *out_table = TRANSLATION_TABLE_5;
                 *out_table_rc = TRANSLATION_TABLE_5_RC;
                 return true;
        case 6: *out_table = TRANSLATION_TABLE_6;
                 *out_table_rc = TRANSLATION_TABLE_6_RC;
                 return true;
        case 9: *out_table = TRANSLATION_TABLE_9;
                 *out_table_rc = TRANSLATION_TABLE_9_RC;
                 return true;
        case 10: *out_table = TRANSLATION_TABLE_10;
                 *out_table_rc = TRANSLATION_TABLE_10_RC;
                 return true;
        case 11: *out_table = TRANSLATION_TABLE_11;
                 *out_table_rc = TRANSLATION_TABLE_11_RC;
                 return true;
        case 12: *out_table = TRANSLATION_TABLE_12;
                 *out_table_rc = TRANSLATION_TABLE_12_RC;
                 return true;
        case 13: *out_table = TRANSLATION_TABLE_13;
                 *out_table_rc = TRANSLATION_TABLE_13_RC;
                 return true;
        case 14: *out_table = TRANSLATION_TABLE_14;
                 *out_table_rc = TRANSLATION_TABLE_14_RC;
                 return true;
        case 15: *out_table = TRANSLATION_TABLE_15;
                 *out_table_rc = TRANSLATION_TABLE_15_RC;
                 return true;
        case 16: *out_table = TRANSLATION_TABLE_16;
                 *out_table_rc = TRANSLATION_TABLE_16_RC;
                 return true;
        case 21: *out_table = TRANSLATION_TABLE_21;
                 *out_table_rc = TRANSLATION_TABLE_21_RC;
                 return true;
        case 22: *out_table = TRANSLATION_TABLE_22;
                 *out_table_rc = TRANSLATION_TABLE_22_RC;
                 return true;
        case 23: *out_table = TRANSLATION_TABLE_23;
                 *out_table_rc = TRANSLATION_TABLE_23_RC;
                 return true;
        case 24: *out_table = TRANSLATION_TABLE_24;
                 *out_table_rc = TRANSLATION_TABLE_24_RC;
                 return true;
        case 25: *out_table = TRANSLATION_TABLE_25;
                 *out_table_rc = TRANSLATION_TABLE_25_RC;
                 return true;
        case 26: *out_table = TRANSLATION_TABLE_26;
                 *out_table_rc = TRANSLATION_TABLE_26_RC;
                 return true;
    }
    return false;
}

#endif
