#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"

#define N_STO 3

typedef struct array {
    double * matrix;
    int row;
    int col;
} Array;

typedef struct array_two {
    double * matrix;
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

    arr->matrix = mat;
    arr->row = x;
    arr->col = y;

    return arr;
}

Array_t * c_arr_doub(int x, int y, int w, int z) 
{   
    Array_t * arr = malloc(sizeof(Array_t));

    double * mat = malloc(sizeof(double) * (x * y * w * z));
    memset(mat, 0, sizeof(double) * (x * y * w * z));

    arr->matrix = mat;
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
            (*S)->matrix[i * (*S)->row + j] = buf[0]; // what if di and dj are > than 1 ???

            cint1e_kin_cart(buf, shls, atm, natm, bas, nbas, env);
            (*T)->matrix[i * (*T)->row + j] = buf[0];

            cint1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env);
            (*V)->matrix[i * (*V)->row + j] = buf[0];
            free(buf);

            (*H)->matrix[i * (*T)->row + j] = (*T)->matrix[i * (*T)->row + j] + (*V)->matrix[i * (*V)->row + j];

            for (int k = 0; k < nbas; k++) {
                for (int l = 0; l < nbas; l++) {
                    shls[2] = k; dk = CINTcgto_cart(k, bas);
                    shls[3] = l; dl = CINTcgto_cart(l, bas);

                    buf = malloc(sizeof(double) * di * dj * dk * dl);
                    cint2e_cart(buf, shls, atm, natm, bas, nbas, env, NULL);
                    // (*two)->matrix[(i * (*two)->x) + (j * (*two)->y) + (k * (*two)->w)]
                    // how ??
                    free(buf);
                }
            }
        }
    }
}

void find_X(double ** S, double *** X, double *** X_dag) 
{
    // sval, U = np.linalg.eig(S)
    // # X = np.matmul(U,np.linalg.inv(s**0.5)) # canonical orthonogalization
    // X = np.matmul(U, np.diag(sval ** (-0.5)))  # symmetric orthogonalization
    // Xdag = X.T
    // return X, Xdag
    return NULL;
}

void calc_F(int nbas, double ** P, double **** two, double ** H, double *** G, double *** F) 
{
    *G = c_arr(nbas, nbas);
    *F = c_arr(nbas, nbas);

    for (int mu = 0; mu < nbas; mu++) {
        for (int nu = 0; nu < nbas; nu++) {
            for (int la = 0; la < nbas; la++) {
                for (int sig = 0; sig < nbas; sig++) {
                    (*G)[mu][nu] += P[la][sig] * (two[mu][nu][sig][la] - 0.5 * two[mu][la][sig][nu]);
                }
            }

            (*F)[mu][nu] += (*G)[mu][nu] + (*H)[mu][nu];
        }
    }
}

void calc_Fprime(double ** F, double ** X, double ** X_dag, double *** Fprime) 
{
    return NULL;
}

void diag_F(double ** Fprime, double ** X, double *** C, double *** epsilon)
{
    // U = np.linalg.eig(Fprime)[1]
    // Udag = np.transpose(U)
    // f = np.matmul(np.matmul(Udag, Fprime), U)
    // Cprime = U
    // epsilon = f
    // C = np.matmul(X, Cprime)
    // return C, epsilon
    return NULL;
}

double ** calc_P(int nbas, int nelec, double ** C) 
{
    double ** P = c_arr(nbas, nbas);
    for (int mu = 0; mu < nbas; mu++) {
        for (int nu = 0; nu < nbas; nu++) {
            for (int i = 0; i < (nelec / 2); i++) {
                P[mu][nu] += 2.0 * C[mu][i] * C[nu][i];
            }
        }
    }
    return P;
}

double f_delta(int nbas, double ** P, double ** P_old) 
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

double RHF(int natm, int nbas, int * atm, int * bas, double * env, int imax, double conv)
{   
    double ** S, ** T, ** V, ** H;
    double **** two;
    integrals(natm, nbas, atm, bas, env, &S, &T, &V, &H, &two);

    double ** X, ** X_dag;
    find_X(S, &X, &X_dag);

    // P is identity
    double ** P = c_arr(nbas, nbas);
    for (int i = 0; i < nbas; i++) {
        for (int j = 0; j < nbas; j++) {
            if (i == j) {
                P[i][j] = 1.0;
            }
        }
    }

    int nelec = 8;

    int i = 0;
    double delta = 1.0;

    double ** P, ** Pold;
    while (delta > conv && i < imax) {
        Pold = P;
        
        double ** G, ** F;
        calc_F(nbas, P, two, H, &G, &F);

        double ** Fprime;
        calc_Fprime(F, X, X_dag, &Fprime);

        double ** C, ** epsilon;
        diag_F(Fprime, X, &C, &epsilon);

        P = calc_P(nbas, nelec, C);

        delta = f_delta(nbas, P, Pold);
        i++;
    }

    if (delta > conv && i == imax) {
        printf("did not converge\n");
    }

    double E0 = 0.0;
    for (int mu = 0; mu < nbas; mu++) {
        for (int nu = 0; nu < nbas; nu++) {
            E0 += 0.5 * P[mu][nu] * (H[mu][nu] + F[mu][nu]);
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

	int * atm = malloc(sizeof(int) * natm * ATM_SLOTS);
	int * bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
	double * env = malloc(sizeof(double) * 10000);

    FILE * file = fopen("/u/jpmedina/libcint/molecules/h2/basis.txt", "r");
 	read_arrays(file, natm, nbas, &atm, &bas, &env);


}
