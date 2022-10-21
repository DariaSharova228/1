#include "matrix.h"
#define EPS 1e-14
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
void get_block_b(int i, int m, int k, int l, double *b, double *block) {
    int rows;
    rows = (i < k ? m : l);
    for (int i1 = 0; i1 < rows; i1++) {
    	block[i1] = b[i * m + i1];
    } 
}
void put_block_b(int i, int m, int k, int l, double *b, double *block) {
    int rows;
    rows = (i < k ? m : l);
    for (int i1 = 0; i1 < rows; i1++) {
    	b[i * m + i] = block[i1];
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
void print_block(double* A, int n1, int n2) {
    for(int i = 0; i < n1; i++){
        for(int j = 0; j < n2; j++) {
            printf("%10.3e", A[i*n2 + j]);
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
    double sum = 0.;
    for (int i = 0; i < n; i++) {
        sum += fabs(X[i]);
    }
    return sum;
}
double norm_block(double* A, int m) {// блок m*m
    double max = 0.;
    double sum = 0.;
    for(int j = 0; j < m; j++) {
        sum = 0.;
        for(int i = 0; i < m; i++) {
            sum += fabs(A[i*m + j]);
        }
        if(sum > max) {
            max = sum;
        }
    }
    return max;
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
int ind_of_max_block(double matrixnorm, int i, int m, double *block) {
    int indmax = -1;
    for(int j = i; j < m; j++) {
        if (!(fabs(block[j * m + i]) < EPS * matrixnorm)) {
            if (indmax == -1) indmax = j;
            else indmax = (fabs(block[indmax * m + i]) - fabs(block[j * m + i]) < EPS * matrixnorm ? j : indmax);
        }
    }
    return indmax;
}
void change_row_block(int i1, int i2, int m, double *block) {
    double buf = 0;
    for(int j = 0; j < m; j++) {
        buf = block[i1 * m + j];
        block[i1 * m + j] = block[i2 * m + j];
        block[i2 * m + j] = buf;
    }
}
void copy(double* block1, double* block2, int m) {
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            block1[i * m + j] = block2[i * m + j];
        }
    }
}
int ind_of_min_matrix(int j, double* A, double normmatrix, double* block1, double* block2, double* block3, int n, int m, int k, int l){
    double curr_norm;
    double min_norm = 0.;
    int min_norm_i = -1;
    for(int i = j; i < k; i++) {//цикл для выбора главного элементa
        get_block(i, j, n, m, k, l, A, block1);
        if(inverse_block(normmatrix, block1, block2, m) == -1){
            return -1;
        }
        if(i == j) {
            min_norm = norm_block(block2, m);
            min_norm_i = i;
        }
        else{
            curr_norm = norm_block(block2, m);
            if(curr_norm < min_norm) {
                min_norm = curr_norm;
                min_norm_i = i;
                copy(block3, block2, m);
                //копировать блок2 в блок3
            }
        }
    }
    return min_norm_i;

}
void change_row_matrix(int i1, int i2, int n, int m, int k, int l, double* A, double* block1, double* block2, double* B){
    for(int j1 = i1; j1 < k; j1++) {
            get_block(i1, j1, n, m, k, l, A, block1);
            get_block(i2, j1, n, m, k, l, A, block2);
            put_block(i2, j1, n, m, k, l, A, block1);
            put_block(i1, j1, n, m, k, l, A, block2);
    }
    get_block(i1, k, n, m, k, l, A, block1);
    get_block(i2, k, n, m, k, l, A, block2);
    put_block(i1, k, n, m, k, l, A, block2);
    put_block(i2, k, n, m, k, l, A, block1);
    get_block_b(i1, m, k, l, B, block1);
    get_block_b(i2, m, k, l, B, block2);
    put_block_b(i2, m, k, l, B, block1); 
    put_block_b(i1, m, k, l, B, block2);
}
/*int inverse_block(double normmatrix, double* block, double* inv_block, int m) {
    int indmax;
    double buf = 0;
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            if(i == j) {
                inv_block[i * m + j] = 1;
            }
            else {
                inv_block[i * m + j] = 0;
            }
        }
    }
    for(int j = 0; j < m; j++) {
        indmax = ind_of_max_block(normmatrix, j, m, block);
        if(indmax != j) {
            change_row_block(j, indmax, m, block);
            change_row_block(j, indmax, m, inv_block);
        }
        if (fabs(block[j * m + j]) < EPS * normmatrix) return -1;
        for(int i = j; i < m; i++) {
            block[j * m + i] = block[j * m + i] / block[j * m + j];
        }
        for(int i = j; i < m; i++) {
            inv_block[j * m + i] = inv_block[j * m + i] / block[j * m + j];
        }
        for(int i = 0; i < m; i++) {
            if (i != j) {
                buf = block[i * m + j];
                for(int j1 = j; j1 < m; j1++) {
                    block[i * m + j1] -= block[j * m + j1] * buf;
                }
                for(int j1 = 0; j1 < m; j1++) {
                    inv_block[i * m + j1] -= inv_block[j * m + j1] * buf;
                }
            }
        }
    }
    return 0;

}*/
int inverse_block(double matrixnorm, double *block, double *block1, int m) {
    int k = 0, i = 0, j = 0;
    double div = 0, buf = 0;
    int indmax = 0;

    for(i = 0; i < m * m; i++) block1[i] = 0;
    for(i = 0; i < m; i++) block1[i * m + i] = 1;

    for(k = 0; k < m; k++) {
        indmax = ind_of_max_block(matrixnorm, k, m, block);
        if (indmax < k) return -1;
        if (indmax > k) {
            change_row_block(k, indmax, m, block);
            change_row_block(k, indmax, m, block1);
        }

        if (fabs(block[k * m + k]) < EPS * matrixnorm) return -1;
        div = 1 / block[k * m + k];

        for(i = k; i < m; i++) 
            block[k * m + i] *= div;
        for(i = 0; i < m; i++)
            block1[k * m + i] *= div;

        for(i = 0; i < m; i++) {
            if (i != k) {
                buf = block[i * m + k];
                for(j = k; j < m; j++) {
                    block[i * m + j] -= block[k * m + j] * buf;
                }
                for(j = 0; j < m; j++) {
                    block1[i * m + j] -= block1[k * m + j] * buf;
                }
            }
        }
    }

    return 1;
}

void mult_blocks(double *block1, double *block2, double *resblock, int n, int m1, int m2) { //block1{m1*m} * block2{m*m2}
    int i = 0, j = 0, k = 0, s = 0, r2 = 0, r3 = 0, l2 = 0, l3 = 0;
    double s00 = 0, s01 = 0, s02 = 0, s10 = 0, s11 = 0, s12 = 0, s20 = 0, s21 = 0, s22 = 0;
    double *b1 = nullptr, *b2 = nullptr, *b3 = nullptr;
    r2 = m1 / 3;
    r3 = m2 / 3;
    l2 = m1 - 3 * r2;
    l3 = m2 - 3 * r3;
    for(i = 0; i < m1 * m2; i++) resblock[i] = 0;
    for(i = 0; i < l2; i++) {
        b1 = block1 + i * n;
        for(j = 0; j < l3; j++)                 
            for(k = 0; k < n; k++)
                resblock[i * m2 + j] += block1[i * n + k] * block2[k * m2 + j];
        for(j = 0; j < r3; j++) {              
            s00 = 0;
            s01 = 0;
            s02 = 0;
            b2 = block2 + l3 + j * 3;
            b3 = resblock + l3 + i * m2 + j * 3;
            for(k = 0; k < n; k++) {
                s00 += b1[k] * b2[k * m2];
                s01 += b1[k] * b2[k * m2 + 1];
                s02 += b1[k] * b2[k * m2 + 2];
            }
            b3[0] += s00;
            b3[1] += s01;
            b3[2] += s02;
        }
    }

    for(j = 0; j < l3; j++) {                   
        b2 = block2 + j;
        for(k = 0; k < r2; k++) {
            s00 = 0;
            s10 = 0;
            s20 = 0;
            b1 = block1 + (l2 + k * 3) * n;
            b3 = resblock + (l2 + k * 3) * m2 + j;
            for(s = 0; s < n; s++) {
                s00 += b1[s] * b2[s * m2];
                s10 += b1[n + s] * b2[s * m2];
                s20 += b1[2 * n + s] * b2[s * m2];
            }
            b3[0] += s00;
            b3[m2] += s10;
            b3[2 * m2] += s20;
        }
    }
    
    for(j = 0; j < r2; j++) {                     
        for(k = 0; k < r3; k++) {
            b1 = block1 + l2 * n + j * 3 * n;
            b2 = block2 + l3 + k * 3;
            b3 = resblock + l2 * m2 + l3 + j * 3 * m2 + k * 3;
            s00 = 0;
            s01 = 0;
            s02 = 0;
            s10 = 0;
            s11 = 0;
            s12 = 0;
            s20 = 0;
            s21 = 0;
            s22 = 0;
            for(s = 0; s < n; s++) {
                s00 += b1[s] * b2[s * m2]; 
                s01 += b1[s] * b2[s * m2 + 1];
                s02 += b1[s] * b2[s * m2 + 2];
                s10 += b1[n + s] * b2[s * m2];
                s11 += b1[n + s] * b2[s * m2 + 1];
                s12 += b1[n + s] * b2[s * m2 + 2];
                s20 += b1[2 * n + s] * b2[s * m2];
                s21 += b1[2 * n + s] * b2[s * m2 + 1];
                s22 += b1[2 * n + s] * b2[s * m2 + 2];
            }
            
            b3[0] += s00;
            b3[1] += s01;
            b3[2] += s02;
            b3[m2] += s10;
            b3[m2 + 1] += s11;
            b3[m2 + 2] += s12;
            b3[2 * m2] += s20;
            b3[2 * m2 + 1] += s21;
            b3[2 * m2 + 2] += s22;
        }
    }
}
void substract_block(double *block1, double *block2, int n, int m) {
    int i = 0, j = 0;
    for(i = 0; i < n; i++)
        for(j = 0; j < m; j++)
            block1[i * m + j] -= block2[i * m + j];    
}
int solution(double* A, double* B, double *X, int n, int m, double* block1, double* block2, double* inv_block, double* block3){
    int k = n/m;
    int l = n - m * k;
    //double curr_norm;
    int min_norm_i;
    double normmatrix = norm_block(A, n);
    for(int p = 0; p < k; p++){//шаг
    min_norm_i = ind_of_min_matrix(p, A, normmatrix, block1, block2, inv_block, n, m, k, l);
    if (min_norm_i == -1) return -1;
    change_row_matrix(p, min_norm_i, n, m, k, l, A, block1, block2, B);
    for(int i1 = p; i1 < k; i1++) {                                
            get_block(i1, p, n, m, k, l, A, block1);        
            mult_blocks(inv_block, block1, block2, m, m, m);
            put_block(i1, p, n, m, k, l, A, block2);
    }
        if (l != 0) { //m x l
            get_block(p, k, n, m, k, l, A, block1);         
            mult_blocks(inv_block, block1, block2, m, m, l);
            put_block(p, k, n, m, k, l, A, block2);
            
        }
        get_block_b(p, m, k, l, B, block1);
        mult_blocks(inv_block, block1, block2, m, m, 1);
        put_block_b(p, m, k, l, B, block2);
        for(int i = p + 1; i < k; i++) {
            get_block(i, p, n, m, k, l, A, block1);
            for(int j = p; j < k; j++) {
                get_block(p, j, n, m, k, l, A, block2);
                mult_blocks(block1, block2, block3, m, m, m);
                get_block(i, j, n, m, k, l, A, block2);
                substract_block(block2, block3, m, m);
                put_block(i, j, n, m, k, l, A, block2);
            }
            if(l != 0) { //m*l
                get_block(p, k, n, m, k, l, A, block2);
                mult_blocks(block1, block2, block3, m, m, l);
                get_block(i, k, n, m, k, l, A, block2);
                substract_block(block2, block3, m, l);
                put_block(i, k, n, m, k, l, A, block2);
            }
            get_block_b(p, m, k, l, B, block2);
            mult_blocks(block1, block2, block3, m, m, 1);
            get_block_b(i, m, k, l, B, block2);
            substract_block(block2, block3, m, 1);
            put_block_b(i, m, k, l, B, block2);
        }
        if (l != 0) { //l x m
            get_block(k, p, n, m, k, l, A, block1);  // l*m      
            for(int j = p; j < k; j++){
                get_block(p, j, n, m, k, l, A, block2);    
                mult_blocks(block1, block2, block3, m, l, m);
                get_block(k, j, n, m, k, l, A, block2);     
                substract_block(block2, block3, l, m);
                put_block(k, j, n, m, k, l, A, block2);
            }
            get_block(p, k, n, m, k, l, A, block2);        
            mult_blocks(block1, block2, block3, m, l, l);
            get_block(k, k, n, m, k, l, A, block2);        
            substract_block(block2, block3, l, l);
            put_block(k, k, n, m, k, l, A, block2);

            get_block_b(p, m, k, l, B, block2);
            mult_blocks(block1, block2, block3, m, l, 1);
            get_block_b(k, m, k, l, B, block2);
            substract_block(block2, block3, l, 1);
            put_block_b(k, m, k, l, B, block2);
        }
    }
    if(l != 0) {
        get_block(k, k, n, m, k, l, A, block1);
        if (inverse_block(normmatrix, block1, inv_block, l) == -1) return -1;
        get_block(k, k, n, m, k, l, A, block1); 
        mult_blocks(inv_block, block1, block3, l, l, l);
        put_block(k, k, n, m, k, l, A, block3);

        get_block_b(k, m, k, l, B, block1);
        mult_blocks(inv_block, block1, block2, l, l, 1);
        put_block_b(k, m, k, l, B, block2);
    }
    for(int k1 = n - 1; k1 >= 0; k1--) {
        X[k1] = B[k1];
        for(int p1 = 0; p1 < k1; p1++) {
            B[p1] = B[p1] - A[p1 * n + k1] * X[k1];
        }
    }
    return 1;

}

