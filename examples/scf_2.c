#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"

#define N_STO 3

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
                    (*two)->m[(int) (i*pow(nbas, 3) + j*pow(nbas, 2) + k*nbas + l)];
                    free(buf);
                }
            }
        }
    }
}

void find_X(Array S, Array * X, Array * X_dag) 
{
    // sval, U = np.linalg.eig(S)
    // # X = np.matmul(U,np.linalg.inv(s**0.5)) # canonical orthonogalization
    // X = np.matmul(U, np.diag(sval ** (-0.5)))  # symmetric orthogonalization
    // Xdag = X.T
    // return X, Xdag
    return NULL;
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

    // free intermediate !!
}

void diag_F(Array Fprime, Array X, Array ** C, Array ** epsilon)
{
    Array * U = c_arr(Fprime.row, Fprime.col);
    // U = np.linalg.eig(Fprime)[1]
    // (how do I get [1] from eig) ??
    // use _CINTdiagonalize ??

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
    memcpy(Cprime->m, U->m, sizeof(U->m));

    *epsilon = c_arr(f->row, f->col);
    // epsilon = f
    memcpy((*epsilon)->m, f->m, sizeof(f->m));

    *C = c_arr(X.row, Cprime->col);
    // C = np.matmul(X, Cprime)
    CINTdgemm_NN(X.row, Cprime->row, Cprime->col, X.m, Cprime->m, (*C)->m);
    
    // return C, epsilon

    // free U, Udag, intermediate, f, Cprime !! 
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
            delta += 0.0; // np.power(P[mu][nu] - P_old[mu][nu], 2.0)
        }
    }
    delta = pow(delta, 0.5) / 2.0;
    return delta;
}

double RHF(int natm, int nbas, int nelec, int * atm, int * bas, double * env, int imax, double conv)
{   
    Array S, T, V, H;
    Array_t two;
    integrals(natm, nbas, atm, bas, env, &S, &T, &V, &H, &two);

    Array X, X_dag;
    find_X(S, &X, &X_dag);

    // P is identity
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
    Array G, F;
    Array Fprime;
    Array C, epsilon;
    while (delta > conv && i < imax) {
        Pold = P;
        
        calc_F(nbas, *P, two, H, &G, &F);
        calc_Fprime(F, X, X_dag, &Fprime);
        diag_F(Fprime, X, &C, &epsilon);

        P = calc_P(nbas, nelec, C);

        delta = f_delta(nbas, *P, *Pold);
        i++;
    }

    if (delta > conv && i == imax) {
        printf("did not converge\n");
    }

    double E0 = 0.0;
    for (int mu = 0; mu < nbas; mu++) {
        for (int nu = 0; nu < nbas; nu++) {
            E0 += 0.5 * P->m[mu * P->row + nu] * (H.m[mu * H.row + nu] + F.m[mu * F.row + nu]);
        }
    }

    // get nuclear energy 
    // is R in env ?

    return E0;
}

int main()
{
    int natm = 2;
	int nbas = 2;
    int nelec = 1; // ?
    int nshells = 1; // ?

	int * atm = malloc(sizeof(int) * natm * ATM_SLOTS);
	int * bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
	double * env = malloc(sizeof(double) * 10000);

    FILE * file = fopen("/u/jpmedina/libcint/molecules/h2/basis.txt", "r");
 	read_arrays(file, natm, nbas, &atm, &bas, &env);


}
