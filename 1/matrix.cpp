#include "matrix.h"
#define EPS 1e-16
void reduce_sum(int p, int *a, int n) {
    static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
    static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
    static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
    static int t_in = 0;
    static int t_out = 0;
    static int *r = nullptr;
    int i = 0;

    if (p <= 1) return ;

    pthread_mutex_lock(&m);
    if (a == nullptr) 
        r = a;
    
    for(i = 0; i < n; i++) 
        r[i] += a[i];

    t_in++;
    if (t_in >= p) {
        t_out = 0;
        pthread_cond_broadcast(&c_in);
    }

    else
        while (t_in < p)
            pthread_cond_wait(&c_in, &m);

    if (r != a) 
        for(i = 0; i < n; i++) a[i] = r[i];
    
    t_out++;
    if(t_out >= p) {
        t_in = 0;
        r = nullptr;
        pthread_cond_broadcast(&c_out);
    }

    else 
        while (t_out < p)
            pthread_cond_wait(&c_out, &m);
    
    pthread_mutex_unlock(&m);
}
int errorHandling(int u, int p, Args* a) {
    int err = 0;
    int i = 0;
    Results *r = a -> res;
    for(i = 0; i < p; i++) {
        
        if(r[i].err < 0) {
            err = r[i].err;
            if(u == 0) {
                switch(r[i].err) {
                    case -1:
                        break;
                    case -2: 
                        printf("Could not open file %s!\n", a[i].name);
                        break;
                    case -3:
                        printf("Could not read file %s!\n", a[i].name);
                        break;
                    case -4:
                        printf("Not enough memory %s!\n", a[i].name);
                        break;
                    default:
                        printf("Unknown error %d %s!\n", r[i].err, a[i].name);
                        break;
                }
            }
        }
    }
    return err;
}

double get_full_time() {
    struct timeval buf;
    gettimeofday(&buf, 0);
    return buf.tv_sec + buf.tv_usec / 1.e6;
}
void* thread_func(void *ptr) {
    Args *a = (Args*) ptr;
    double *A = a -> A;
    double *B = a -> B;
    double *X = a -> X;
    Results *res = a -> res;
    Global_results *g_r = a -> global_res;
    int p = a -> p;
    int u = a -> u;
    int k = a -> k;
    int n = a -> n;
    int l = a -> l;
    int m = a -> m;
    int r = a -> r;
    int s = a -> s;
    double matrixnorm = a -> matrixnorm;
    int indmin = -1;
    int err = 0;
    double tmp = 0.;
    double norm = 0, minnorm = 0, normb = 0;
    double t_cpu = 0, t_tot = 0;
    char *name = a -> name;
    FILE *inp;
    int i = 0, j = 0, h = 0;
    int task = 9;
    double r1 = -1, r2 = -1, t1 = 0, t2 = 0;
    double* block1;
    double* block2;
    double* block3;
    double* inv_block;
    double *inv_block1 = nullptr;
    cpu_set_t cpu;
    CPU_ZERO(&cpu);
    int n_cpus = get_nprocs();
    int cpu_id = n_cpus - 1 - (u % n_cpus);
    CPU_SET(cpu_id, &cpu);
    pthread_t tid = pthread_self();
    pthread_setaffinity_np(tid, sizeof(cpu), &cpu);

    for(i = u * m; i < n; i += p * m) { 
        h = (i + m <= n ? m : n - i);
        memset(A + i * n, 0, h * n * sizeof(double));
        memset(X + i, 0, h * sizeof(double));
    }
    reduce_sum(p);
    if (name) {
        if (u == 0) {
            inp = fopen(name, "r");
            if (!inp) {
                res[u].err = -2;
            }
            else {
                if (read_matrix(inp, A, n) == -1) {
                    res[u].err = -3;
                }
            }
            fclose(inp);
        }
    }
    else {
        fill_matrix(A, n, s, m, u, p);
    }
    reduce_sum(p);
    err = errorHandling(u, p, a);
    if(err != 0) {
        if(u == 0) {
            g_r -> err = err;
        }
        return nullptr;
    }
    initialization_B(A, B, n, m, u, p);
    reduce_sum(p);
    initialization_X(X, n, m, u, p);
    reduce_sum(p);

    if (u == 0) {
        printf("A = \n");
        print_matrix(A, n, m, r);
        printf("\n");
    }
    
    /*block1 = new double[m * m]; 
    inv_block = new double[m * m];
    block2 = new double[m * m]; 
    block3 = new double[m * m];
    
    if ((!block3) || (!block2) || (!block1) || (!inv_block)) {
        if (block1) delete []block1;
        if (inv_block) delete []inv_block;
        if (block2) delete []block2;
        if (block3) delete []block3;
        res[u].err = -4;
        return nullptr;
    }*/
    block1 = new double[m*m];
    if(!block1) {
        res[u].err = -4;
        return nullptr;
    }
    block2 = new double[m*m];
    if(!block2) {
        delete []block1;
        res[u].err = -4;
        return nullptr;
    }
    block3 = new double[m*m];
    if(!block3) {
        delete []block1;
        delete []block2;
        res[u].err = -4;
        return nullptr;
    }
    inv_block = new double[m*m];
    if(!inv_block) {
        delete []block1;
        delete []block2;
        delete []block3;
        res[u].err = -4;
        return nullptr;
    }

    reduce_sum(p);
    err = errorHandling(u, p, a);
    if(err != 0) {
        if(u == 0) {
            g_r -> err = err;
        }
        return nullptr;
    }

    for(i = 0; i < n; i++) {
        norm_matrix(A, &res[u].norm, i, u, p, n);
        reduce_sum(p);
        for(j = 0; j < p; j++) 
            norm += res[j].norm;
        reduce_sum(p);
        if (fabs(norm) > fabs(matrixnorm))
            matrixnorm = fabs(norm);
        res[u].norm = 0;
        norm = 0;
    }
    //printf("%d        %lf\n", u, matrixnorm);

    if (n == 1 && A[0] < EPS) {
        res[u].err = -1;
    }
    reduce_sum(p);
    err = errorHandling(u, p, a);
    if(err != 0) {
        if(u == 0) {
            g_r -> err = err;
        }
        return nullptr;
    }

    t_cpu = get_full_time();
    t_tot = get_full_time();
    for(j = 0; j < k; j++) {
        res[u].ind = ind_of_min_matrix(matrixnorm, j, A, block1, block2, inv_block, &res[u].minnorm, n, m, k, l, u, p);
        res[u].invblock = inv_block;
        reduce_sum(p);
        for(i = 0; i < p; i++) {
            if ((fabs(res[i].minnorm) - fabs(minnorm) < EPS * matrixnorm || indmin == -1) && res[i].ind != -1) {
                minnorm = res[i].minnorm;
                indmin = res[i].ind;
                inv_block1 = res[i].invblock;
                //for(j = 0; j < m * m; j++)
                //    invblock[j] = res[i].invblock[j];
            }
        }        
                
        if (indmin == -1) {
            res[u].err = -1;
            delete [] block1;
            delete [] block2;
            delete [] block3;
            delete [] inv_block;
            break;
        }
        change_row_matrix(j, indmin, n, m, k, l, p, u, A, B, block1, block2);
        reduce_sum(p);
        if(!(solution(matrixnorm, A, B, X, block1, block2, inv_block1, block3, n, m, k, l, u, p, j))) {
            res[u].err = -1;
            delete [] block1;
            delete [] block2;
            delete [] block3;
            delete [] inv_block;
            break;
        }
        minnorm = 0;
        indmin = -1;
    }
    reduce_sum(p);
    err = errorHandling(u, p, a);
    if(err != 0) {
        if(u == 0) {
            g_r -> err = err;
        }
        return nullptr;
    }
    if (l > 0) {
        if (!(solution(matrixnorm, A, B, X, block1, block2, inv_block1, block3, n, m, k, l, u, p, j))) {
            res[u].err = -1;
            delete [] block1;
            delete [] block2;
            delete [] block3;
            delete [] inv_block;
        }
        reduce_sum(p);
        err = errorHandling(u, p, a);
        if(err != 0) {
            if(u == 0) {
                g_r -> err = err;
            }
            return nullptr;
        }
    }
    a -> t_cpu = get_full_time() - t_cpu;
    a -> t_tot = get_full_time() - t_tot;    
    if (l == 0) {
        if(u == 0) {
            for (int i = n - 1; i >= 0; --i){
                tmp = B[i];
                for (int j = i + 1; j < n; ++j)
                    tmp -= A[i * n + j] * X[j];
                X[i] = tmp;
            }
        }
        reduce_sum(p);
    }
    /*if (u == 0) {
        printf("A = \n");
        printMatrix(r, n, n, matrix);
        printf("\n");
    }*/
    if (name) {
        if (u == 0) {
            inp = fopen(name, "r");
            if (!inp) {
                res[u].err = -2;
            }
            else {
                if (read_matrix(inp, A, n) == -1) {
                    res[u].err = -3;
                }
            }
            fclose(inp);
        }
    }
    else {
        fill_matrix(A, n, s, m, u, p); 
    }
    reduce_sum(p);
    err = errorHandling(u, p, a);
    if(err != 0) {
        if(u == 0) {
            g_r -> err = err;
        }
        return nullptr;
    }
    initialization_B(A, B, n, m, u, p);
    reduce_sum(p);
    t1 = a -> t_tot;
    t2 = get_full_time();
    normB(B, &normb, n, u, p);
    reduce_sum(p);
    discrepancy1(A, X, B, &r1, n, m, u, p);
    reduce_sum(p);
    r1 = (normb < EPS * matrixnorm ? -1 : (r1 + 1.) / normb);
    discrepancy2(X, &r2, n, u, p);
    reduce_sum(p);
    t2 = get_full_time() - t2;
    if (u == 0) {
        printf("x = \n");
        print_matrix(X, n, 1, r);
        printf("\n");
        printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", "./a.out", task, r1, r2, t1, t2, s, n, m, p);
    }
    delete []block1;
    delete []block2;
    delete []block3;
    delete []inv_block;
    return nullptr;
}
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
    	b[i * m + i1] = block[i1];
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
        case 4: return ((1./((double)(i + j + 1))));
        default: return 0.;
    }
}
void fill_matrix(double* A, int n, int s, int m, int u, int p) {
    int i = 0, j = 0, i1 = 0, h = 0;
    for(i = u * m; i < n; i += p * m) {
        h = (i + m < n ? m : n - i);
        for(i1 = i; i1 < i + h; i1++)
            for(j = 0; j < n; j++) {
                A[i1 * n + j] = 0;
                A[i1 * n + j] = filling_formulas(s, n, i1, j);
            }
    }
}
/*void fill_matrix(double* A, int n, int s){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++)
        {
            A[i * n + j] = 0;
            A[i*n + j] = filling_formulas(s, n, i, j);
        }
    }
}*/
/*void initialization_B(double* B, double* A,int n){
    for(int i = 0; i < n; i++) {
        B[i] = 0. ;
        for(int k = 0; k < n; k += 2){
            B[i] += A[i*n + k]; 
        }
    }
}*/
void initialization_B(double *A, double *B, int n, int m, int u, int p) {
    int i = 0, j = 0, i1 = 0, k = 0, h = 0;
    k = (n + 1) / 2;
    for(i = u * m; i < n; i += p * m) {
        h = (i + m < n ? m : n - i);
        for(i1 = i; i1 < i + h; i1++) {
            B[i1] = 0;
            for(j = 0; j < k; j++)
                B[i1] += A[i1 * n + 2 * j];//возможно неправильно
        }
    }
}
void initialization_X(double *X, int n, int m, int u, int p) {
    int i = 0, i1 = 0, h = 0;
    for(i = u * m; i < n; i += p * m) {
        h = (i + m < n ? m : n - i);
        for(i1 = i; i1 < i + h; i1++) {
            X[i1] = 0;
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
/*double norm_matrix(double* A, int m) {
    double max = 0.;
    double sum = 0.;
    for(int i = 0; i < m; i++) {
        sum = 0.;
        for(int j = 0; j < m; j++) {
            sum += fabs(A[i*m + j]);
        }
        if(sum > max) {
            max = sum;
        }
    }
    return max;
}*/
void norm_matrix(double *matrix, double *matrixnorm, int j, int u, int p, int n) {
    int i = 0;
    double norm = 0;

    for(i = u; i < n; i += p) {
        norm += fabs(matrix[j * n + i]);
    }

    *matrixnorm += norm;
}
double norm_block(double* block, double normmatrix, int m) {
    double norm = 0;
    double sum = 0;
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < m; j++) {
            sum += fabs(block[j * m + i]);
        }
        norm = (fabs(norm) - fabs(sum) < EPS * normmatrix ? sum : norm);
    }
    return norm;
}
/*double discrepancy1(double* A, double* B, double* X, int n) {
    
    int i = 0, j = 0;
    double res = 0, bnorm = 0, s0 = 0, s1 = 0, s2 = 0;
    double matrixnorm = norm_matrix(A, n);
    double *m = nullptr;

    int n3 = n - n % 3;
    for(i = 0; i < n3; i += 3) {
        bnorm += fabs(B[i]) + fabs(B[i + 1]) + fabs(B[i + 2]);
        
        s0 = 0;        
        s1 = 0;
        s2 = 0;
        m = A + i * n;
        
        for(j = 0; j < n3; j += 3) {
            s0 += m[j] * X[j] + m[j + 1] * X[j + 1] + m[j + 2] * X[j + 2];
            s1 += m[n + j] * X[j] + m[n + j + 1] * X[j + 1] + m[n + j + 2] * X[j + 2];
            s2 += m[2 * n + j] * X[j] + m[2 * n + j + 1] * X[j + 1] + m[2 * n + j + 2] * X[j + 2];
        }
        for(j = n3; j < n; j++) {
            s0 += m[j] * X[j];
            s1 += m[n + j] * X[j];
            s2 += m[2 * n + j] * X[j];
        }
        
        res += fabs(s0 - B[i]) + fabs(s1 - B[i + 1]) + fabs(s2 - B[i + 2]);

    }
    for(i = n3; i < n; i++) {
        bnorm += fabs(B[i]);

        s0 = 0; 
        m = A + i * n;
        for(j = 0; j < n3; j += 3) 
            s0 += m[j] * X[j] + m[j + 1] * X[j + 1] + m[j + 2] * X[j + 2];
        
        for(j = n3; j < n; j++) 
            s0 += m[j] * X[j];
        res += fabs(s0 - B[i]);
    }
    
    return (bnorm < EPS * matrixnorm ? -1 : res / bnorm);
}*/
void discrepancy1(double *matrix, double *x, double *b, double *res0, int n, int m, int u, int p) {
    int i = 0, j = 0, h = 0, i1 = 0;
    double res = 0, s = 0;

    for(i = u * m; i < n; i += p * m) {
        h = (i + m < n ? m : i - n);
        for(i1 = i; i1 < i + h; i1++) {
            s = 0;
            for(j = 0; j < n; j++) 
                s += matrix[i1 * n + j] * x[j];
            res += fabs(s - b[i1]);
        }
    }
    
    *res0 += res;
}
/*double discrepancy2(double* X, int n) {
    double sum = 0.;
    for(int i = 0; i < n; i++) {
        sum += fabs(X[i] - ((i + 1) & 1U));
    }
    return sum;
}*/
void discrepancy2(double *x, double *res, int n, int u, int p) {
    int i = 0;
    double res0 = 1.;
    for(i = u; i < n; i += p)
        res0 += fabs(x[i] - ((i + 1) & 1U));
    *res += res0;
}
int ind_of_max_block(double matrixnorm, int i, int m, double *block) {
    int indmax = -1;
    for(int j = i; j < m; j++) {
        if (!(fabs(block[j * m + i]) < EPS * matrixnorm )) {
            if (indmax == -1) indmax = j;
            else indmax = (fabs(block[indmax * m + i]) - fabs(block[j * m + i]) < EPS * matrixnorm? j : indmax);
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
int ind_of_min_matrix(double normmatrix, int j, double* A,double* block1, double* block2, double* block3, double *min, int n, int m, int k, int l, int u, int p){
    double curr_norm;
    double min_norm = 0.;
    int min_norm_i = -1;
    for(int i = j + u; i < k; i += p) {//цикл для выбора главного элементa
        get_block(i, j, n, m, k, l, A, block1);
        if(inverse_block(normmatrix, block1, block2, m) != -1){
            curr_norm = norm_block(block2,normmatrix, m);
            if ((fabs(curr_norm) - fabs(min_norm) < EPS * normmatrix) || (min_norm_i == -1)) {
                    min_norm = curr_norm;
                    min_norm_i = i;
                    copy(block3, block2, m);
            }
        }
    }
    *min = min_norm;
    return min_norm_i;

}
void change_row_matrix(int i1, int i2, int n, int m, int k, int l, int p, int u, double* A, double* B, double* block1, double* block2){
    for(int j1 = i1 + u; j1 < k + 1; j1 += p) {
        get_block(i1, j1, n, m, k, l, A, block1);
        get_block(i2, j1, n, m, k, l, A, block2);
        put_block(i2, j1, n, m, k, l, A, block1);
        put_block(i1, j1, n, m, k, l, A, block2);
        if(j1 == k && l > 0) {
            get_block(i1, k, n, m, k, l, A, block1);
            get_block(i2, k, n, m, k, l, A, block2);
            put_block(i1, k, n, m, k, l, A, block2);
            put_block(i2, k, n, m, k, l, A, block1);
        }
    }
    if(u == 0) {
        get_block_b(i1, m, k, l, B, block1);
        get_block_b(i2, m, k, l, B, block2);
        put_block_b(i2, m, k, l, B, block1); 
        put_block_b(i1, m, k, l, B, block2);
    }
}
int inverse_block(double matrixnorm, double *block, double *block1, int m) {
    int k = 0, i = 0, j = 0;
    double mult = 0, tmp = 0;
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
        mult = 1 / block[k * m + k];

        for(i = k; i < m; i++) 
            block[k * m + i] *= mult;
        for(i = 0; i < m; i++)
            block1[k * m + i] *= mult;

        for(i = 0; i < m; i++) {
            if (i != k) {
                tmp = block[i * m + k];
                for(j = k; j < m; j++) {
                    block[i * m + j] -= block[k * m + j] * tmp;
                }
                for(j = 0; j < m; j++) {
                    block1[i * m + j] -= block1[k * m + j] * tmp;
                }
            }
        }
    }

    return 1;
}

void mult_blocks(double *block1, double *block2, double *resblock, int n, int m1, int m2) { //block1{m1*n} * block2{n*m2}
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
/*void clean_block(double* block, int m) {
    for(int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            block[i * m + j] = 0.;
        }
    }
}*/
int solution(double normmatrix, double* A, double* B, double *X, double* block1, double* block2, double* inv_block, double* block3, int n, int m, int k, int l, int u, int p, int s) {
    double tmp = 0.;
    if(s < k) {
        for(int i1 = s + u; i1 < k + 1; i1 += p) {
            if (l > 0 && i1 == k) {
                get_block(s, k, n, m, k, l, A, block1);       
                mult_blocks(inv_block, block1, block2, m, m, l);
                put_block(s, k, n, m, k, l, A, block2);
            }
            else {
                get_block(s, i1, n, m, k, l, A, block1); 
                mult_blocks(inv_block, block1, block2, m, m, m);
                put_block(s, i1, n, m, k, l, A, block2);
            }
        }
        if (u == 0) {
            get_block_b(s, m, k, l, B, block1);
            mult_blocks(inv_block, block1, block2, m, m, 1);
            put_block_b(s, m, k, l, B, block2);
        }
        reduce_sum(p);
        for(int i = u + s; i < k + 1; i += p) {
            if (i != s && i != k) {
                get_block(i, s, n, m, k, l, A, block1);        //i_s
                for(int j = s; j < k; j++) {
                        get_block(s, j, n, m, k, l, A, block2);    //s_j
                        mult_blocks(block1, block2, block3, m, m, m);
                        get_block(i, j, n, m, k, l, A, block2);     //i_j
                        substract_block(block2, block3, m, m);
                        put_block(i, j, n, m, k, l, A, block2);
                    } 
                if (l > 0) { //m x l
                    get_block(s, k, n, m, k, l, A, block2);    //s_k+1
                    mult_blocks(block1, block2, block3, m, m, l);
                    get_block(i, k, n, m, k, l, A, block2);     //i_k+1
                    substract_block(block2, block3, m, l);
                    put_block(i, k, n, m, k, l, A, block2);
                }
                get_block_b(s, m, k, l, B, block2);
                mult_blocks(block1, block2, block3, m, m, 1);
                get_block_b(i, m, k, l, B, block2);
                substract_block(block2, block3, m, 1);
                put_block_b(i, m, k, l, B, block2);
            }
            if (l > 0 && i == k) { //l x m
                get_block(k, s, n, m, k, l, A, block1);        //k+1_s
                    for(int j = s; j < k; j++) {
                        get_block(s, j, n, m, k, l, A, block2);    //s_j
                        mult_blocks(block1, block2, block3, m, l, m);
                        get_block(k, j, n, m, k, l, A, block2);     //k+1_j
                        substract_block(block2, block3, l, m);
                        put_block(k, j, n, m, k, l, A, block2);
                    }
                    get_block(s, k, n, m, k, l, A, block2);        //s_k+1
                    mult_blocks(block1, block2, block3, m, l, l);
                    get_block(k, k, n, m, k, l, A, block2);         //k+1_k+1
                    substract_block(block2, block3, l, l);
                    put_block(k, k, n, m, k, l, A, block2);

                    get_block_b(s, m, k, l, B, block2);
                    mult_blocks(block1, block2, block3, m, l, 1);
                    get_block_b(k, m, k, l, B, block2);
                    substract_block(block2, block3, m, 1);
                    put_block_b(k, m, k, l, B, block2);
            }
        }
    }
    if(l > 0 && s == k) {
        get_block(k, k, n, m, k, l, A, block1);
        if (inverse_block(normmatrix, block1, inv_block, l) == -1) return -1;
        if(u == 0) {
            get_block(k, k, n, m, k, l, A, block1);
            mult_blocks(inv_block, block1, block3, l, l, l);
            put_block(k, k, n, m, k, l, A, block3);

            get_block_b(k, m, k, l, B, block1);
            mult_blocks(inv_block, block1, block2, l, l, 1);
            put_block_b(k, m, k, l, B, block2);
        }
    }
    if(u == 0) {
    for (int i = n - 1; i >= 0; --i)
	{
		tmp = B[i];
		for (int j = i + 1; j < n; ++j)
			tmp -= A[i * n + j] * X[j];
		X[i] = tmp;
	}
    }
    
    reduce_sum(p);
    return 1;
}
void normB(double *b, double *norm, int n, int u, int p) {
    int i = 0;
    double res = 0;
    for(i = u; i < n; i += p)
        res += fabs(b[i]);
    *norm += res;
}
/*int solution(double* A, double* B, double *X, double* block1, double* block2, double* inv_block, double* block3, int n, int m, int k, int l, int u, int p, int s){
    int p;
    double tmp = 0.;
    //double s = 0.;
    //double curr_norm;
    int min_norm_i;
    double normmatrix = norm_matrix(A, n);
    for( p = 0; p < k; p++){//шаг
    min_norm_i = ind_of_min_matrix(normmatrix, p, A, block1, block2, inv_block, n, m, k, l);
    if (min_norm_i == -1) return -1;
    change_row_matrix(p, min_norm_i, n, m, k, l, A, block1, block2, B);
    for(int i1 = p; i1 < k; i1++) {                                
            get_block(p, i1, n, m, k, l, A, block1);    
            mult_blocks(inv_block, block1, block2, m, m, m);
            put_block(p, i1, n, m, k, l, A, block2);
            //clean_block(block1, m);
            //clean_block(block2, m);
    }
        if (l != 0) { //m x l
            get_block(p, k, n, m, k, l, A, block1);   
            mult_blocks(inv_block, block1, block2, m, m, l);
            put_block(p, k, n, m, k, l, A, block2);
            //clean_block(block1, m);
            //clean_block(block2, m);
            
        }
        get_block_b(p, m, k, l, B, block1);
        mult_blocks(inv_block, block1, block2, m, m, 1);
        put_block_b(p, m, k, l, B, block2);
        //clean_block(block1, m);
        //clean_block(block2, m);
        for(int i = p + 1; i < k; i++) {
            get_block(i, p, n, m, k, l, A, block1);
            for(int j = p; j < k; j++) {
                get_block(p, j, n, m, k, l, A, block2);
                mult_blocks(block1, block2, block3, m, m, m);
                get_block(i, j, n, m, k, l, A, block2);
                substract_block(block2, block3, m, m);
                put_block(i, j, n, m, k, l, A, block2);
                //clean_block(block2, m);
                //clean_block(block3, m);
            }
            if(l != 0) { //m*l
                get_block(p, k, n, m, k, l, A, block2);
                mult_blocks(block1, block2, block3, m, m, l);
                get_block(i, k, n, m, k, l, A, block2);
                substract_block(block2, block3, m, l);
                put_block(i, k, n, m, k, l, A, block2);
                //clean_block(block2, m);
                //clean_block(block3, m);
            }
            get_block_b(p, m, k, l, B, block2);
            mult_blocks(block1, block2, block3, m, m, 1);
            //clean_block(block2, m);
            get_block_b(i, m, k, l, B, block2);
            substract_block(block2, block3, m, 1);
            put_block_b(i, m, k, l, B, block2);
            //clean_block(block2, m);
            //clean_block(block3, m);
        }
        //clean_block(block1, m);
        if (l > 0) { //l x m
            get_block(k, p, n, m, k, l, A, block1);  // l*m      
            for(int j = p; j < k; j++){
                get_block(p, j, n, m, k, l, A, block2);
                mult_blocks(block1, block2, block3, m, l, m);
                get_block(k, j, n, m, k, l, A, block2);     
                substract_block(block2, block3, l, m);
                put_block(k, j, n, m, k, l, A, block2);
                //clean_block(block2, m);
                //clean_block(block3, m);
            }
            get_block(p, k, n, m, k, l, A, block2);
            mult_blocks(block1, block2, block3, m, l, l);
            get_block(k, k, n, m, k, l, A, block2);        
            substract_block(block2, block3, l, l);
            put_block(k, k, n, m, k, l, A, block2);
           //clean_block(block2, m);
            //clean_block(block3, m);
            get_block_b(p, m, k, l, B, block2);
            mult_blocks(block1, block2, block3, m, l, 1);
            get_block_b(k, m, k, l, B, block2);
            substract_block(block2, block3, l, 1);
            put_block_b(k, m, k, l, B, block2);
            //clean_block(block2, m);
            //clean_block(block3, m);
            //clean_block(block1, m);
        }
    }
    if(l > 0) {
        get_block(k, k, n, m, k, l, A, block1);
        if (inverse_block(normmatrix, block1, inv_block, l) == -1) return -1;
        get_block(k, k, n, m, k, l, A, block1);
        mult_blocks(inv_block, block1, block3, l, l, l);
        put_block(k, k, n, m, k, l, A, block3);

        get_block_b(k, m, k, l, B, block1);
        mult_blocks(inv_block, block1, block2, l, l, 1);
        put_block_b(k, m, k, l, B, block2);
    }
    for (int i = n - 1; i >= 0; --i)
	{
		tmp = B[i];
		for (int j = i + 1; j < n; ++j)
			tmp -= A[i * n + j] * X[j];
		X[i] = tmp;
	}
    return 1;
}*/
