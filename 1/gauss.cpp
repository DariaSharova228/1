#include "matrix.h"

int main(int argc, char *argv[]) {
    Args *args = nullptr;
    Results *res = nullptr;
    Global_results g_r;
    int n = 0, m = 0, r = 0, s = 0, p = 0, u = 0, l = 0, k = 0;
    double *A = nullptr;
    double *X = nullptr;
    double *B = nullptr;
    char *name = nullptr;
    //double r1 = -1, r2 = -1, t1 = 0, t2 = 0;
    //FILE *inp;

    if (!((argc == 6 || argc == 7) 
    && sscanf(argv[1], "%d", &n)
    && sscanf(argv[2], "%d", &m)
    && sscanf(argv[3], "%d", &p)
    && sscanf(argv[4], "%d", &r)
    && sscanf(argv[5], "%d", &s)
    && (s >= 0) && (s <= 4)) || (argc == 6 && s == 0)) {
        printf("%s Usage: n m p r s [file]\n", argv[0]);
        return 0;
    }
    
    A = new double[n * n];
    if (!A) {
        printf("MEMORY ERROR!!!\n");
        return 0;
    }
    X = new double[n];
    if (!X) {
        printf("MEMORY ERROR!!!\n");
        delete []A;
        return 0;
    }

    B = new double[n];
    if (!B) {
        printf("MEMORY ERROR!!!\n");
        delete []A;
        delete []X;
        return 0;
    }

    if (s == 0) {
        name = argv[6];
    }

    res = new Results[p];
    if (res == nullptr) {
        printf("MEMORY ERROR!!!\n");
        delete []A;
        delete []X;
        delete []B;
        return 0;
    }

    args = new Args[p];
    if (args == nullptr) {
        printf("MEMORY ERROR!!!\n");
        delete []A;
        delete []X;
        delete []B;
        delete []res;
        return 0;
    }

    k = n / m;
    l = n % m;
    for(u = 0; u < p; u++) {
        args[u].A = A;
        args[u].B = B;
        args[u].X = X;
        args[u].name = name;
        args[u].argv_0 = argv[0];
        args[u].n = n;
        args[u].m = m;
        args[u].r = r;
        args[u].s = s;
        args[u].u = u;
        args[u].p = p;
        args[u].k = k;
        args[u].l = l;
        args[u].res = res;
        args[u].global_res = &g_r;
    }

    for(u = 1; u < p; u++) {
        if(pthread_create(&args[u].t_id, 0, thread_func, args + u)) {
            printf("COULD NOT CREATE PTHREAD %d!\n", u);
            delete []A;
            delete []X;
            delete []B;
            delete []res;
            delete []args;
            abort();
        }
    }

    thread_func(args + 0);
    for(u = 1; u < p; u++) {
        pthread_join(args[u].t_id, 0);
    }

    if (g_r.err < 0) printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d P = %d\n", argv[0], 9, -1., -1., 0., 0., s, n, m, p);

    delete []A;
    delete []X;
    delete []B;
    delete []args;
    delete []res;
    return 0;
}
