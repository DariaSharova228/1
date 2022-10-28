#include <stdio.h>
void put_block(int i, int j, int n, int m, int k, int l, double* A, double* block);
void get_block(int i, int j, int n, int m, int k, int l, double* A, double* block);
void get_block_b(int i, int m, int k, int l, double *b, double *block);
int read_matrix(FILE* inp, double* A, int n);
void print_block(double* A, int n1, int n2);
void print_matrix(double *A, int n, int m, int r);
double filling_formulas(int s, int i, int j, int n);
void fill_matrix(double* A, int n, int s);
void initialization_B(double* B, double* A,int n);
void multiplication_block(double* A, double* B, double* C, int rows1, int cols1, /*int rows2,*/ int cols2);
double norm_vector(double* X, int n);
double norm_matrix(double* A, int m);  
double norm_block(double* block, double normmatrix, int m);
double discrepancy1(double* A, double* B, double* X, int n);
double discrepancy2(double* X, int n);
int ind_of_max_block(double matrixnorm, int j, int m, double *block);
void change_row_block(int i1, int i2, int m, double *block);
void copy(double* block1, double* block2, int m);
int inverse_block(double matrixnorm, double* block, double* inv_block, int m);
int ind_of_min_matrix(double normmatrix, int j, double* A, double* block1, double* block2, double* block3, int n, int m, int k, int l);
void change_row_matrix(int i1, int i2, int n, int m, int k, int l, double* A, double* block1, double* block2, double* B);
//void mult_blocks(double *a, double *b, double *res, int m2, int m1, int m3, int m);
void mult_blocks(double *block1, double *block2, double *resblock, int n, int m1, int m2);//block1{m1*n} * block2{n*m2}
int solution(double* A, double* B, double *X, int n, int m, double* block1, double* block2, double* inv_block, double* block3);
