#include <stdio.h>
#include <time.h>
#include "matrix.h"
int main(int argc, char *argv[]){
    double* A;
    double* B;
    double* X;
    double* block1;
    double* block2;
    double* block3;
    double* inv_block;
    int res = 0;
    int n = 0, m = 0, r = 0, s = 0;
    char* name_file;
    int task = 9;
    double r1 = 0, r2 = 0, t1 = 0, t2 = 0;
    FILE* inp;
    if(!((argc != 5) || (argc != 6))){
        return 0;
    }
    if(!(sscanf(argv[1], "%d", &n) == 1 && 
    sscanf(argv[2], "%d", &m) == 1 && 
    sscanf(argv[3], "%d", &r) == 1 && 
    ((sscanf(argv[4], "%d", &s) == 1 &&
    (s >= 0 && s <= 4)) || (s == 0 && argc == 5))
    )) {
        return 0;
    }
    A = new double[n*n];
    if(!A){
        return 0;
    }
    B = new double[n];
    if(!B)
    {
        delete []A;
        return 0;
    }
    X = new double[n];
    if(!X){
        delete []A;
        delete []B;
        return 0;
    }
    if(s == 0) {
        name_file = argv[5];
        inp = fopen(name_file, "r");
        if(!inp){
            delete []A;
            delete []B;
            delete []X;
            return 0;
        }
        if(read_matrix(inp, A, n) == 0 ) {
            printf("ошибка в файле с матрицей\n");\
            fclose(inp);
            delete []A;
            delete []B;
            delete []X;
            return 0;
        }
        initialization_B(B, A, n);
        fclose(inp);
    }
    else {
        fill_matrix(A, n, s);
        initialization_B(B, A, n);
    }
    printf("A =\n");
  print_matrix(A, n, n, r);
  printf("\n");
  printf("B =\n");
  print_matrix(B, n, 1, r);
  printf("\n");
  //int k = n/m;
  //int l = n - m * k;
  block1 = new double[m*m];
  block2 = new double[m*m];
  block3 = new double[m*m];
  inv_block = new double [m*m];
  for(int i = 0; i < n; i++) {
    X[i] = 0;
  }
  /*double normmatrix = norm_block(A, n);
  get_block(0, 0, n, m, k, l, A, block1);
  print_block(block1, m, m);
  printf("\n");
  get_block_b(0, m, k, l , B, block2);
  print_block(block2, m, 1);
  printf("\n");
  mult_blocks(block1, block2, block3, m, m, 1);
  print_block(block3, m, 1);
  printf("\n");*/
  /*get_block(0, 0, n, m, k, l, A, block1);
  print_block(block1, m, m);
  printf("\n");
  get_block(1, 0, n, m, k, l, A, block2);
  print_block(block2, l, m);
  printf("\n");
  mult_blocks(block2, block1, block3, m, l, m);
  print_block(block3, l, m);
  printf("\n");*/
  /*printf("B =\n");
  print_matrix(B, n, 1, r);*/
    t1 = clock();
    res = solution(A, B, X, n, m, block1, block2, inv_block, block3);
    //printf("A =\n");
  //print_matrix(A, n, n, r);
  //printf("\n");
  //printf("B =\n");
  //print_matrix(B, n, 1, r);
  //printf("\n");
    t1 = (clock() - t1)/ CLOCKS_PER_SEC;
    if (res != 1){
        r1 = -1;
        r2 = -1;
        printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], task, r1, r2, t1, t2, s, n, m);
        delete []A;
        delete []B;
        delete []X;
        delete [] block1;
        delete [] block2;
        delete [] block3;
        delete [] inv_block;
        return 0;
    }
    printf("X =\n");
    print_matrix(X, n, 1, r);
    if(s == 0) {
        name_file = argv[5];
        inp = fopen(name_file, "r");
        if(!inp){
            delete []A;
            delete []B;
            delete []X;
            delete [] block1;
            delete [] block2;
            delete [] block3;
            delete [] inv_block;
            return 0;
        }
        if(read_matrix(inp, A, n) == 0 ) {
            printf("ошибка в файле с матрицей\n");\
            fclose(inp);
            delete []A;
            delete []B;
            delete []X;
            delete [] block1;
            delete [] block2;
            delete [] block3;
            delete [] inv_block;
            return 0;
        }
        initialization_B(B, A, n);
        fclose(inp);
    }
    else {
        fill_matrix(A, n, s);
        initialization_B(B, A, n);
    }
    t2 = clock();
    r1 = discrepancy1(A, B, X, n);
    r2 = discrepancy2(X, n);
    t2 = (clock() - t2)/CLOCKS_PER_SEC;
    printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], task, r1, r2, t1, t2, s, n, m);
    delete []A;
    delete []B;
    delete []X;
    delete []block1;
    delete []block2;
    delete []block3;
    delete []inv_block;
    return 0;
}
