#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"

#include "f2c.h"
#include "blaswrap.h"

#define N_STO 3

int __enzyme_autodiff(void *, ...);

int cint1e_ovlp_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_kin_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_nuc_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

extern int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a, 
    integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *info);

extern int enzyme_dup;
extern int enzyme_out;
extern int enzyme_const;

typedef struct array {
    double * m;
    int row;
    int col;
} Array;

typedef struct Array_two {
    double * m;
    int x;
    int y;
    int w;
    int z;
} Array_two;

void free_arr(Array * arr) {
    free(arr->m);
    free(arr);
}

void free_arr_two(Array_two * arr) {
    free(arr->m);
    free(arr);
}

void read_arrays(FILE * file, int natm, int nbas, int ** atm, int ** bas, double ** env)
{
	double num;
	int number;
	int count = 0;

	while (count < natm * ATM_SLOTS) {
        fscanf(file, "%d", &number);
		(*atm)[count] = number;
		count++;
	}

	count = 0;
	while (count < nbas * BAS_SLOTS) {
        fscanf(file, "%d", &number);
		(*bas)[count] = number;
		count++;
	}

	count = 0;
	while (fscanf(file, "%lf", &num) != EOF) {
		(*env)[count] = num;
		count++;
	}
	fclose(file);
}

Array * c_arr(int x, int y) 
{   
    Array * arr = malloc(sizeof(Array));

    double * mat = malloc(sizeof(double) * (x * y));
    memset(mat, 0, sizeof(double) * (x * y));

    arr->m = mat;
    arr->row = x;
    arr->col = y;

    return arr;
}

Array_two * c_arr_doub(int x, int y, int w, int z) 
{   
    Array_two * arr = malloc(sizeof(Array_two));

    double * mat = malloc(sizeof(double) * (x * y * w * z));
    memset(mat, 0, sizeof(double) * (x * y * w * z));

    arr->m = mat;
    arr->x = x;
    arr->y = y;
    arr->w = w;
    arr->z = z;

    return arr;
}

void dcopy(Array * dest, Array source) {
    for (int i = 0; i < source.col * source.row; i++) {
        dest->m[i] = source.m[i];
    }
}

void matmult(int n, double * A, double * B, double * C) {
    // cij = aik bkj
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

void transpose(int n, double * C, double * Ct) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Ct[i * n + j] = C[j * n + i];
        }
    }
}

void integrals(int natm, int nbas, int * atm, int * bas, double * env, Array ** S, Array ** T, Array ** V, Array ** H, Array_two ** two) 
{   
    *S = c_arr(nbas, nbas);
    *T = c_arr(nbas, nbas);
    *V = c_arr(nbas, nbas);
    *H = c_arr(nbas, nbas);
    *two = c_arr_doub(nbas, nbas, nbas, nbas);

    int di, dj, dk, dl;
    int shls[4];
    
    double * buf;
    for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
            shls[0] = i; di = CINTcgto_cart(i, bas);
            shls[1] = j; dj = CINTcgto_cart(j, bas);

            buf = malloc(sizeof(double) * di * dj);
            cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);
            (*S)->m[i * (*S)->row + j] = buf[0];

            cint1e_kin_cart(buf, shls, atm, natm, bas, nbas, env);
            (*T)->m[i * (*T)->row + j] = buf[0];

            cint1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env);
            (*V)->m[i * (*V)->row + j] = buf[0];
            free(buf);

            (*H)->m[i * (*T)->row + j] = (*T)->m[i * (*T)->row + j] + (*V)->m[i * (*V)->row + j];

            for (int k = 0; k < nbas; k++) {
                for (int l = 0; l < nbas; l++) {
                    shls[2] = k; dk = CINTcgto_cart(k, bas);
                    shls[3] = l; dl = CINTcgto_cart(l, bas);

                    buf = malloc(sizeof(double) * di * dj * dk * dl);
                    cint2e_cart(buf, shls, atm, natm, bas, nbas, env, NULL);
                    (*two)->m[(int) (i*pow(nbas, 3) + j*pow(nbas, 2) + k*nbas + l)] = buf[0];
                    free(buf);
                }
            }
        }
    }
}

void find_X(Array S, Array ** X, Array ** X_dag) 
{
    Array * eig = c_arr(1, S.col);
    Array * U = c_arr(S.row, S.col);
    dcopy(U, S);

    char jobz = 'V', uplo = 'U';
    integer lda = 2, n = 2, info = 8, lwork = 6;
    double w[2], work[6];

    dsyev_(&jobz, &uplo, &n, U->m, &lda, w, work, &lwork, &info);

    if(info > 0) {
        printf("The algorithm failed to compute eigenvalues.\n");
        exit(1);
    }

    eig->m = w;

    Array * diag_eig = c_arr(S.col, S.col);
    for (int i = 0; i < S.col; i++) {
        for (int j = 0; j < S.col; j++) {
            if (i == j) {
                diag_eig->m[2*i + j] = pow(eig->m[i], -0.5);
            }
        }
    }

    *X = c_arr(U->row, diag_eig->row);
    matmult(U->row, U->m, diag_eig->m, (*X)->m);

    *X_dag = c_arr((*X)->row, (*X)->col);
    transpose((*X)->row, (*X)->m, (*X_dag)->m);
}

void calc_F(int nbas, Array P, Array_two two, Array H, Array ** G, Array ** F) 
{
    *G = c_arr(nbas, nbas);
    *F = c_arr(nbas, nbas);

    for (int mu = 0; mu < nbas; mu++) {
        for (int nu = 0; nu < nbas; nu++) {
            for (int la = 0; la < nbas; la++) {
                for (int sig = 0; sig < nbas; sig++) {
                    (*G)->m[mu * (*G)->row + nu] += P.m[la * P.row + sig] 
                        * (two.m[(int) (mu*pow(nbas, 3) + nu*pow(nbas, 2) + sig*nbas + la)] 
                        - 0.5 * two.m[(int) (mu*pow(nbas, 3) + la*pow(nbas, 2) + sig*nbas + nu)]);
                }
            }

            (*F)->m[mu * (*F)->row + nu] += (*G)->m[mu * (*G)->row + nu] + H.m[mu * H.row + nu];
        }
    }
}

void calc_Fprime(Array F, Array X, Array X_dag, Array ** Fprime) 
{
    Array * inter = c_arr(X_dag.row, F.col);
    matmult(X_dag.row, X_dag.m, F.m, inter->m);

    *Fprime = c_arr(inter->row, X.col);
    matmult(inter->row, inter->m, X.m, (*Fprime)->m);

    free_arr(inter);
}

void diag_F(Array Fprime, Array X, Array ** C, Array ** epsilon)
{
    Array * U = c_arr(Fprime.row, Fprime.col);
    dcopy(U, Fprime);

    char jobz = 'V', uplo = 'U';
    integer lda = 2, n = 2, info = 8, lwork = 6;
    double w[2], work[6];

    dsyev_(&jobz, &uplo, &n, U->m, &lda, w, work, &lwork, &info);

    if(info > 0) {
        printf("The algorithm failed to compute eigenvalues.\n");
        exit(1);
    }

    Array * Udag = c_arr(U->row, U->col);
    transpose(Udag->row, U->m, Udag->m);
    
    Array * inter = c_arr(Udag->row, Fprime.col);
    matmult(Udag->row, Udag->m, Fprime.m, inter->m);
    
    Array * f = c_arr(inter->row, U->col);
    matmult(inter->row, inter->m, U->m, f->m);

    Array * Cprime = c_arr(U->row, U->col);
    dcopy(Cprime, *U);

    *epsilon = c_arr(f->row, f->col);
    dcopy(*epsilon, *f);

    *C = c_arr(X.row, Cprime->col);
    matmult(X.row, X.m, Cprime->m, (*C)->m);
    
    free_arr(U);
    free_arr(Udag);
    free_arr(inter);
    free_arr(f);
    free_arr(Cprime);
}

Array * calc_P(int nbas, int nelec, Array C) 
{
    Array * P = c_arr(nbas, nbas);
    for (int mu = 0; mu < nbas; mu++) {
        for (int nu = 0; nu < nbas; nu++) {
            for (int i = 0; i < (nelec / 2); i++) {
                // printf("%lf\n", 2.0 * C.m[mu * C.row + i] * C.m[nu * C.row + i]);
                P->m[mu * P->row + nu] += 2.0 * C.m[mu * C.row + i] * C.m[nu * C.row + i];
            }
        }
    }

    return P;
}

double f_delta(int nbas, Array P, Array P_old) 
{
    double delta = 0;
    for (int mu = 0; mu < nbas; mu++) {
        for (int nu = 0; nu < nbas; nu++) {
            delta += pow(P.m[mu * P.row + nu] - P_old.m[mu * P.row + nu], 2.0);
        }
    }
    delta = pow(delta, 0.5) / 2.0;

    return delta;
}

double norm(int * atm, double * env, int i, int j) 
{
    double xi = env[atm[i*6 + 1]];
    double xj = env[atm[j*6 + 1]];

    double yi = env[atm[i*6 + 1] + 1];
    double yj = env[atm[j*6 + 1] + 1];

    double zi = env[atm[i*6 + 1] + 2];
    double zj = env[atm[j*6 + 1] + 2];

    return pow(pow(xi - xj, 2) + pow(yi - yj, 2) + pow(zi - zj, 2), 0.5); 
}

Array * RHF(int natm, int nbas, int nelec, int * atm, int * bas, double * env, int imax, double conv)
{   
    Array * S, * T, * V, * H;
    Array_two *two;
    integrals(natm, nbas, atm, bas, env, &S, &T, &V, &H, &two);

    Array * X, * X_dag;
    find_X(*S, &X, &X_dag);

    Array * P = c_arr(nbas, nbas);
    for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
            if (i == j) {
                P->m[i * P->row + j] = 1.0;
            }
        }
    }

    int i = 0;
    double delta = 1.0;

    Array * Pold;
    Array * G, * F;
    Array * Fprime;
    Array * C, * epsilon;
    
    while (delta > conv && i < imax) {
        Pold = P;
        
        calc_F(nbas, *P, *two, *H, &G, &F);
        calc_Fprime(*F, *X, *X_dag, &Fprime);
        diag_F(*Fprime, *X, &C, &epsilon);

        P = calc_P(nbas, nelec, *C);

        delta = f_delta(nbas, *P, *Pold);
        i++;
    }

    if (delta > conv && i == imax) {
        printf("did not converge\n");
        exit(1);
    }

    return P;
}

void energy(double * E, int natm, int nbas, int * atm, int * bas, double * env, Array * P)
{
    Array * S, * T, * V, * H;
    Array_two *two;
    integrals(natm, nbas, atm, bas, env, &S, &T, &V, &H, &two);

    Array * G, * F;    
    calc_F(nbas, *P, *two, *H, &G, &F);

    double E0 = 0.0;
    for (int mu = 0; mu < nbas; mu++) {
        for (int nu = 0; nu < nbas; nu++) {
            E0 += 0.5 * P->m[mu * P->row + nu] * (H->m[mu * H->row + nu] + F->m[mu * F->row + nu]);
        }
    }

    double Enuc = 0.0;
    for (int i = 0; i < natm; i++) {
        for (int j = 0; j < natm; j++) {
            if (i > j) {
                Enuc += ((double) (atm[i*6 + 0] * atm[j*6 + 0]))/(norm(atm, env, i, j));
            }
        }
    }

    *E = Enuc + E0;
}

void energy_diff(double * E, double * dE, int natm, int nbas, int * atm, int * bas, double * env, double * denv, Array * P) {
	__enzyme_autodiff((void *) energy,
        enzyme_dup, E, dE,	
        enzyme_const, natm,
        enzyme_const, nbas,
        enzyme_const, atm,
        enzyme_const, bas, 
        enzyme_dup, env, denv,
        enzyme_const, P);
}

int main()
{
    int natm = 2;
    int nbas = 2;
    
    int nelec = 2; // ?
    int nshells = 1; // ?

    int * atm = malloc(sizeof(int) * natm * ATM_SLOTS);
    int * bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
    double * env = malloc(sizeof(double) * 10000);

    FILE * file = fopen("/u/jpmedina/libcint/molecules/h2/basis.txt", "r");
    read_arrays(file, natm, nbas, &atm, &bas, &env);

    // RHF
    int imax = 20;
    double conv = 0.000001;

    Array * P = RHF(natm, nbas, nelec, atm, bas, env, imax, conv);

    double E;
    double dE = 1.0;
    double * denv = malloc(sizeof(double) * 10000);

    energy_diff(&E, &dE, natm, nbas, atm, bas, env, denv, P);
    printf("E: %lf\n", E);
    printf("denv:\n");
    for (int k = 28; k < 34; k++) {
        printf("%f ", denv[k]);
    }
    printf("\n");

    printf("finite diff:\n");
    
    double grad;
    double E1, E2;
    double h = 0.000001;
    
    for (int k = 28; k < 34; k++) {
        env[k] += h;
        energy(&E1, natm, nbas, atm, bas, env, P);
        env[k] -= 2.0*h;
        energy(&E2, natm, nbas, atm, bas, env, P);
        env[k] += h;

        grad = (E1 - E2)/(2.0*h);

        printf("%lf ", grad);
    }
    printf("\n");

    return 0;
}
