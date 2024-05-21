#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"

#define N_STO 3

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

double ** c_arr(int x, int y) 
{
    double ** arr = malloc(sizeof(double *) * x);

    for (int i = 0; i < x; i++) {
        arr[i] = malloc(sizeof(double) * y);
        memset(arr[i], 0, sizeof(double) * y);
    }

    return arr;
}

double **** c_arr_doub(int x, int y, int w, int z) 
{
    double **** arr = malloc(sizeof(double ***) * x);

    for (int i = 0; i < x; i++) {
        arr[i] = malloc(sizeof(double **) * y);
        for (int j = 0; j < y; j++) {
            arr[i][j] = malloc(sizeof(double *) * w);
            for (int k = 0; k < w; k++) {
                arr[i][j][k] = malloc(sizeof(double) * z);
                memset(arr[i][j][k], 0, sizeof(double) * z);
            }
        }
    }

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


void integrals(int natm, int nbas, int * atm, int * bas, double * env, double *** S, double *** T, double *** V, double *** H, double ***** two) 
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
            (*S)[i][j] = buf[0]; // what if di and dj are > than 1 ???

            cint1e_kin_cart(buf, shls, atm, natm, bas, nbas, env);
            (*T)[i][j] = buf[0];

            cint1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env);
            (*V)[i][j] = buf[0];
            free(buf);

            (*H)[i][j] = (*T)[i][j] + (*V)[i][j];

            for (int k = 0; k < nbas; k++) {
                for (int l = 0; l < nbas; l++) {
                    shls[2] = k; dk = CINTcgto_cart(k, bas);
                    shls[3] = l; dl = CINTcgto_cart(l, bas);

                    buf = malloc(sizeof(double) * di * dj * dk * dl);
                    cint2e_cart(buf, shls, atm, natm, bas, nbas, env, NULL);
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
