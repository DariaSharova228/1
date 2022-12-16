#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <time.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <sys/resource.h>

#define EPS 1e-16

class Results {
    public:
        double norm = 0;
        double minnorm = -1;
        double *invblock = nullptr;
        int ind = -1;
        int err = 0;

        Results() = default;
        ~Results() = default;
};

class Global_results {
    public:
        int err = 0;

        Global_results() = default;
        ~Global_results() = default;
};

class Args {
    public:
        Results *res = nullptr;
        double *A = nullptr;
        double *B = nullptr;
        double *X = nullptr;
        int n = 0;
        int m, l, k, r, s;
        int u; //number of a thread
        int p = 0;
        double matrixnorm = 0;
        char *name = nullptr;
        pthread_t t_id = -1;
        char *argv_0 = nullptr;
        double t_cpu = 0;
        double t_tot = 0;
        Global_results *global_res = nullptr;
        
        Args() = default;
        ~Args() = default;
};
double fabs(double a);
void init_block0(int n, int m, double *block);
void reduce_sum(int p, int *a = nullptr, int n = 0);
double get_full_time();
void normB(double *b, double *norm, int n, int u, int p);
void* thread_func(void *ptr);
void put_block(int i, int j, int n, int m, int k, int l, double* A, double* block);
void get_block(int i, int j, int n, int m, int k, int l, double* A, double* block);
void get_block_b(int i, int m, int k, int l, double *b, double *block);
void put_block_b(int i, int m, int k, int l, double *b, double *block);
int read_matrix(FILE* inp, double* A, int n);
void print_block(double* A, int n1, int n2);
void print_matrix(double *A, int n, int m, int r);
double filling_formulas(int s, int n, int i, int j);
void fill_matrix(double* A, int n, int s, int m, int u, int p);
void initialization_B(double *A, double *B, int n, int m, int u, int p);
void initialization_X(double *X, int n, int m, int u, int p);
void multiplication_block(double* A, double* B, double* C, int rows1, int cols1, /*int rows2,*/ int cols2);
double norm_vector(double* X, int n);
void norm_matrix(double *matrix, double *matrixnorm, int j, int u, int p, int n);
double norm_block(double* block, double normmatrix, int m);
void discrepancy1(double *matrix, double *x, double *b, double *res0, int n, int m, int u, int p);
void discrepancy2(double *x, double *res, int n, int u, int p);
int ind_of_max_block(double matrixnorm, int j, int m, double *block);
void change_row_block(int i1, int i2, int m, double *block);
void copy(double* block1, double* block2, int m);
int inverse_block(double matrixnorm, double* block, double* inv_block, int m);
int ind_of_min_matrix(double normmatrix, int j, double* A,double* block1, double* block2, double* block3, double *min, int n, int m, int k, int l, int u, int p);
void change_row_matrix(int i1, int i2, int n, int m, int k, int l, int p, int u, double* A, double* B, double* block1, double* block2);
//void mult_blocks(double *a, double *b, double *res, int m2, int m1, int m3, int m);
void mult_blocks(double *block1, double *block2, double *resblock, int n, int m1, int m2);//block1{m1*n} * block2{n*m2}
int solution(double normmatrix, double* A, double* B, double* block1, double* block2, double* inv_block, double* block3, int n, int m, int k, int l, int u, int p, int s);
