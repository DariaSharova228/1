#include "matrix.h"
void get_block(int i, int j, int n, int m, int k, int l, double *A, double *block){
    int rows, cols;
    rows = (i == k? l : m);
    cols = (j == k? l : m);
    for(int i1 = 0; i1 < rows; i1++)
    {
        for(int j1 = 0; j1 < cols; j1++)
        {
            block[i1 * cols + j1] = A[n * (i * m + i1) + j * m + j1];
        }
    }
}
void put_block(int i, int j, int n, int m, int k, int l, double *A, double *block){
    int rows, cols;
    rows = (i == k? l : m);
    cols = (j == k? l : m);
    for(int i1 = 0; i1 < rows; i1++){
        for(int j1 = 0; j1 < cols; j1++){
            A[n * (i * m + i1) + j * m + j1] = block[i1 * cols + j1];
        }
    }
}
double fabs(double a) {
    return (a > 0 ? a : -a);
}
int read_matrix(FILE* inp, double* A, int n){
    int i;
    int j;
    double a;
    for(i = 0; i < n; i++) {
        for (j = 0; j < n; j++){
            if(fscanf(inp, "%lf", &a) != 1) {
                return 0;
            }
            A[i*n + j] = a;
        }
    }
    return 1;
}
void print_matrix(double* A, int n, int m, int r){
  int i1 = (n > r ? r : n);
  int j1 = (m > r ? r : m);
  for (int i = 0; i < i1; i++)
    {
      for(int j = 0; j < j1; j++)
        {
          printf("%10.3e", A[i*m + j]);
        }
      printf("\n");
    }
}
double filling_formulas(int s, int n, int i, int j) {
    switch(s) {
        case 1: return ( n - (i > j ? i : j));
        case 2: return (i < j ? i + 1: j + 1);
        case 3: return(i - j < 0 ? j - i : i - j);
        case 4: return ((1/((double)(i + j + 1))));
        default: return 0.;
    }
}
void fill_matrix(double* A, int n, int s){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
        {
            A[i*n + j] = filling_formulas(s, n, i, j);
        }
    }
}
void initialization_B(double* B, double* A,int n){
    for(int i = 0; i < n; i++) {
        B[i] = 0. ;
        for(int k = 0; k < n; k += 2){
            B[i] += A[i*n + k]; 
        }
    }
}
//void multiplication_matr(double* A, double* B)
void multiplication_block(double* A, double* B, double* C, int rows1, int cols1,/* int rows2,*/ int cols2) { //C=AB
    for(int i = 0; i < rows1; i++) {
        for(int j = 0; j < cols2; j++) {
            for(int s = 0; s < cols1; s++) {
                C[i*cols2 + j] += A[i*cols1 + s] * B[s*cols2 + j]; 
            }
        }
    }
}
double norm_vector(double* X, int n){
    int sum = 0.;
    for (int i = 0; i < n; i++) {
        sum += fabs(X[i]);
    }
    return sum;
}
double discrepancy1(double* A, double* B, double* X, int n){
    double sum1 = 0.;
    double sum2 = 0.;
    double norm;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++){
            sum1 += A[i*n +j] * X[j];
        }
        sum2 += fabs(sum1 - B[i]);
        sum1 = 0.;
    }
    norm = norm_vector(B,n);
    return (norm > 0 ? (sum2/norm) : 10000);
}
/*double discrepancy1(double* A, double* B, double* X, int n){
    double sum = 0.;
    int k = n/3;

}*/
double discrepancy2(double* X, int n) {
    double sum = 0.;
    for(int i = 0; i < n; i++) {
        sum += fabs(X[i] - (i % 2));
    }
    return sum;
}
int solution(double *X, int n){
  for(int i = 0; i < n; i+=2){
    X[i] = 1;
  }
  for(int i = 1; i < n; i +=2 ) {
    X[i] = 0;
  }
  return 1;
}

