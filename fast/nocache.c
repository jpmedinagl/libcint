#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "cint.h"

#define ATM_SLOTS 6
#define BAS_SLOTS 8

int int1e_nuc_cart(double *out, int *dims, int *shls, int *atm, int natm,
    int *bas, int nbas, double *env, CINTOpt *opt, double *cache);

int main() {
    int natm = 3;
    int nbas = 12;

    int atm[] = {
        8, 20, 1, 23, 0, 0,
        1, 24, 1, 27, 0, 0,
        1, 28, 1, 31, 0, 0
    };

    int bas[] = {
        0,  0,  5,  1,  0, 32, 37,  0,
        0,  0,  1,  1,  0, 42, 43,  0,
        0,  0,  1,  1,  0, 44, 45,  0,
        0,  1,  3,  1,  0, 46, 49,  0,
        0,  1,  1,  1,  0, 52, 53,  0,
        0,  2,  1,  1,  0, 54, 55,  0,
        1,  0,  3,  1,  0, 56, 59,  0,
        1,  0,  1,  1,  0, 62, 63,  0,
        1,  1,  1,  1,  0, 64, 65,  0,
        2,  0,  3,  1,  0, 56, 59,  0,
        2,  0,  1,  1,  0, 62, 63,  0,
        2,  1,  1,  1,  0, 64, 65,  0
    };

    double env[] = {
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, -0.2104232716416691, 0.0, 0.0, 0.0, 0.841692897594064,
        -1.4797241525927651, 0.0, 0.0, 0.841692897594064, 1.4797241525927651,
        0.0, 2266.1767785, 340.87010191, 77.363135167, 21.47964494,
        6.6589433124, -4.472205830377161, -8.064133436980633, -11.86821744420372,
        -11.804535745439507, -4.680635259503801, 0.80975975668, 2.1566623979211927,
        0.25530772234, 0.907429696013769, 17.721504317, 3.863550544, 1.0480920883,
        6.643459018549882, 5.2670542098593485, 2.293965466701481, 0.27641544411,
        0.5847056269797676, 1.2, 3.590017508385001, 13.010701, 1.9622572,
        0.44453796, 0.5795583105047568, 0.9831856491001469, 1.1193051215867884,
        0.12194962, 0.5213751919783473, 0.8, 2.207226371076266
    };

    int shls[4] = {0, 0, 0, 0};

    CINTOpt *opt = NULL;

    double *buf = NULL;
    int di, dj;

    printf("buf\n");
    clock_t start = clock();
    for (int i = 0; i < nbas; ++i) {
        for (int j = 0; j < nbas; ++j) {
            shls[0] = i;
            shls[1] = j;

            di = CINTcgto_cart(i, bas);
            dj = CINTcgto_cart(j, bas);

            buf = (double *)malloc(sizeof(double) * di * dj);

            int1e_nuc_cart(
                buf,
                NULL,
                shls,
                atm,
                natm,
                bas,
                nbas,
                env,
                opt,
                NULL
            );
        }
    }
    clock_t end = clock();

    double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Elapsed time: %.6f seconds\n", elapsed_time);

    return 0;
}
