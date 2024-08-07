#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"

// #include "f2c.h"
// #include "blaswrap.h"

#define N_STO 3

int cint1e_ovlp_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_kin_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_nuc_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int _CINTdiagonalize(int n, double *diag, double *diag_off1, double *eig, double *vec);

void CINTdgemm_NN(FINT m, FINT n, FINT k, double *a, double *b, double *c);

void CINTdmat_transpose(double *a_t, double *a, FINT m, FINT n);

// extern int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a, 
//     integer *lda, doublereal *w, doublereal *work, integer *lwork, 
// 	integer *info);

typedef struct array {
    double * m;
    int row;
    int col;
} Array;

typedef struct array_two {
    double * m;
    int x;
    int y;
    int w;
    int z;
} Array_t;

void free_arr(Array * arr) {
    free(arr->m);
    free(arr);
}

void free_arr_t(Array_t * arr) {
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

Array_t * c_arr_doub(int x, int y, int w, int z) 
{   
    Array_t * arr = malloc(sizeof(Array_t));

    double * mat = malloc(sizeof(double) * (x * y * w * z));
    memset(mat, 0, sizeof(double) * (x * y * w * z));

    arr->m = mat;
    arr->x = x;
    arr->y = y;
    arr->w = w;
    arr->z = z;

    return arr;
}

// STO_NG ?? information already in env, atm, bas

// void STO_NG(int nbas, double * zeta, double *** alpha, double *** D)
// {
//     *alpha = c_arr(nbas, N_STO);
//     *D = c_arr(nbas, N_STO);
    
//     for (int i = 0; i < nbas, i++) {
//         alpha
//     }
// }


void integrals(int natm, int nbas, int * atm, int * bas, double * env, Array ** S, Array ** T, Array ** V, Array ** H, Array_t ** two) 
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
    // sval, U = np.linalg.eig(S)

    Array * diag = c_arr(S.row, 1);
    Array * diag_off = c_arr(S.row - 1, 1);
    for (int i = 0; i < S.row; i++) {
        for(int j = 0; j < S.col; j++) {
            if (i == j) {
                diag->m[i] = S.m[i * S.row + j];
            } else if (i - 1 == j) {
                diag_off->m[i] = S.m[i * S.row + j];
            }
        }
    }
    
    // ERROR 1
    Array * eig = c_arr(1, S.col);
    Array * U = c_arr(S.row, S.col);
    _CINTdiagonalize(S.row, diag->m, diag_off->m, eig->m, U->m);

    printf("Eig Vec\n");
    for (int i = 0; i < S.row; i++) {
        printf("%lf\n", eig->m[i]);
        for (int j = 0; j < S.row; j++) {
            printf("%lf ", U->m[i * U->row + j]);
        }
        printf("\n");
    }
    printf("\n");

    Array * diag_eig = c_arr(S.col, S.col);
    for (int i = 0; i < S.col; i++) {
        for (int j = 0; j < S.col; j++) {
            if (i == j) {
                diag_eig->m[i] = pow(eig->m[i], -0.5);
            }
        }
    }

    // X = np.matmul(U, np.diag(sval ** (-0.5)))  # symmetric orthogonalization
    *X = c_arr(U->row, diag_eig->row);
    CINTdgemm_NN(U->row, diag_eig->row, diag_eig->col, U->m, diag_eig->m, (*X)->m);

    // Xdag = X.T
    *X_dag = c_arr((*X)->row, (*X)->col);
    CINTdmat_transpose((*X_dag)->m, (*X)->m, (*X)->row, (*X)->col);
    
    // return X, Xdag
}

void calc_F(int nbas, Array P, Array_t two, Array H, Array ** G, Array ** F) 
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
    // Fprime = np.matmul(np.matmul(Xdag, F), X)

    Array * inter = c_arr(X_dag.row, F.col);
    CINTdgemm_NN(X_dag.row, F.row, F.col, X_dag.m, F.m, inter->m);

    *Fprime = c_arr(inter->row, X.col);
    CINTdgemm_NN(inter->row, X.row, X.col, inter->m, X.m, (*Fprime)->m);

    // free intermediate
    free_arr(inter);
}

void diag_F(Array Fprime, Array X, Array ** C, Array ** epsilon)
{
    // Array * U = c_arr(Fprime.row, Fprime.col);
    // U = np.linalg.eig(Fprime)[1]

    // Array * eig = c_arr(1, Fprime.col);
    // don't need the eigenvalues
    Array * diag = c_arr(Fprime.row, 1);
    Array * diag_off = c_arr(Fprime.row - 1, 1);
    for (int i = 0; i < Fprime.row; i++) {
        for(int j = 0; j < Fprime.col; j++) {
            if (i == j) {
                diag->m[i] = Fprime.m[i * Fprime.row + j];
            } else if (i - 1 == j) {
                diag_off->m[i] = Fprime.m[i * Fprime.row + j];
            }
        }
    }

    Array * eig = c_arr(1, Fprime.col);
    Array * U = c_arr(Fprime.row, Fprime.col);
    _CINTdiagonalize(Fprime.row, diag->m, diag_off->m, eig->m, U->m);

    Array * Udag = c_arr(U->row, U->col);
    // Udag = np.transpose(U)
    CINTdmat_transpose(Udag->m, U->m, U->row, U->col);
    
    Array * inter = c_arr(Udag->row, Fprime.col);
    // np.matmul(Udag, Fprime)
    CINTdgemm_NN(Udag->row, Fprime.row, Fprime.col, Udag->m, Fprime.m, inter->m);
    
    Array * f = c_arr(inter->row, U->col);
    // f = np.matmul(np.matmul(Udag, Fprime), U)
    CINTdgemm_NN(inter->row, U->row, U->col, inter->m, U->m, f->m);

    Array * Cprime = c_arr(U->row, U->col);
    // Cprime = U
    memcpy(Cprime->m, U->m, U->row * U->col);

    // ERROR 2: memcpy is not correct

    *epsilon = c_arr(f->row, f->col);
    // epsilon = f
    memcpy((*epsilon)->m, f->m, f->row * f->col);

    *C = c_arr(X.row, Cprime->col);
    // C = np.matmul(X, Cprime)
    CINTdgemm_NN(X.row, Cprime->row, Cprime->col, X.m, Cprime->m, (*C)->m);
    
    // return C, epsilon

    // free U, Udag, intermediate, f, Cprime
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

double RHF(int natm, int nbas, int nelec, int * atm, int * bas, double * env, int imax, double conv)
{   
    Array *S, *T, *V, *H;
    Array_t *two;
    integrals(natm, nbas, atm, bas, env, &S, &T, &V, &H, &two);
    printf("S T V\n");
    for (int i = 0; i < nbas * nbas; i++) {
        printf("%lf %lf %lf\n", S->m[i], T->m[i], V->m[i]);
    }
    printf("\n");

    printf("two");
    for (int i = 0; i < nbas * nbas * nbas * nbas; i++) {
        if (i % 4 == 0) {
            printf("\n");
        }
        printf("%lf ", two->m[i]);
    }
    printf("\n\n");

    Array *X, *X_dag;
    find_X(*S, &X, &X_dag);
    printf("X Xdag\n");
    for (int i = 0; i < nbas * nbas; i++) {
        printf("%lf %lf\n", X->m[i], X_dag->m[i]);
    }
    printf("\n");

    // P is identity
    Array * P = c_arr(nbas, nbas);
    for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
            if (i == j) {
                P->m[i * P->row + j] = 1.0;
            }
        }
    }
    
    // for (int i = 0; i < nbas; i++) {
    //     for (int j = 0; j < nbas; j++) {
    //         printf("%lf ", P->m[i * P->row + j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    int i = 0;
    double delta = 1.0;

    Array * Pold;
    Array *G, *F;
    Array *Fprime;
    Array *C, *epsilon;
    
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
    }

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

    double Etot = Enuc + E0;
    return Etot;
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

    double Etot = RHF(natm, nbas, nelec, atm, bas, env, imax, conv);

    printf("Etot: %lf \n", Etot);

    return 0;
}
