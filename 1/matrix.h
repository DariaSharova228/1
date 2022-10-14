#include <stdio.h>
void put_block(int i, int j, int n, int m, int k, int l, double* A, double* block);
void get_block(int i, int j, int n, int m, int k, int l, double* A, double* block);
int read_matrix(FILE* inp, double* A, int n);
void print_matrix(double *A, int n, int m, int r);
double filling_formulas(int s, int i, int j, int n);
void fill_matrix(double* A, int n, int s);
void initialization_B(double* B, double* A,int n);
void multiplication_block(double* A, double* B, double* C, int rows1, int cols1, /*int rows2,*/ int cols2);
double norm_vector(double* X, int n);
double discrepancy1(double* A, double* B, double* X, int n);
double discrepancy2(double* X, int n);
int solution( double *X, int n);
