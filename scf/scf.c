#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cint.h"

#include "f2c.h"
#include "blaswrap.h"

#define N_STO 3

extern int enzyme_dup;
extern int enzyme_out;
extern int enzyme_const;

int __enzyme_autodiff(void *, ...);

typedef int (* int1e)(double *buf, int *shls, 
					int *atm, int natm, int *bas, int nbas, double *env);

typedef int (* int2e)(double *buf, int *shls, 
					int *atm, int natm, int *bas, int nbas, double *env, CINTOpt *opt);

typedef FINT (* cgto)(const FINT bas_id, const FINT *bas);

int cint1e_ovlp_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_ovlp_sph(double *buf, int *shls, 
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_kin_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_kin_sph(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_nuc_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_nuc_sph(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

extern int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a, 
    integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *info);

double * c_arr(int x, int y) 
{
    double * mat = malloc(sizeof(double) * (x * y));
    memset(mat, 0, sizeof(double) * (x * y));
    return mat;
}

double * c_arr_doub(int x, int y, int w, int z) 
{   
    double * mat = malloc(sizeof(double) * (x * y * w * z));
    memset(mat, 0, sizeof(double) * (x * y * w * z));
    return mat;
}

void print_arr(int n, int size, double * A) {
    for (int i = 0; i < pow(n, size); i++) {
        printf("%lf ", A[i]);
        for (int p = 1; p < size; p++) {
            if ((i + 1) % (int) pow(n, p) == 0) {
                printf("\n");
            }
        }
    }
}

void dcopy(int n, double * dest, double * source) {
    for (int i = 0; i < n * n; i++) {
        dest[i] = source[i];
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

double * integral1e(int natm, int nbas, int nshells, int * atm, int * bas, double * env, int coord, int type) {
    double * R = c_arr(nshells, nshells);

    int1e func;
    cgto intcgto;

    if (coord == 0) {
        // CARTERSIAN
        intcgto = CINTcgto_cart;

        if (type == 0) {
            func = cint1e_ovlp_cart;
        } else if (type == 1) {
            func = cint1e_kin_cart;
        } else if (type == 2) {
            func = cint1e_nuc_cart;
        }
    } else if (coord == 1) {
        // SPHERIC
        intcgto = CINTcgto_spheric;

        if (type == 0) {
            func = cint1e_ovlp_sph;
        } else if (type == 1) {
            func = cint1e_kin_sph;
        } else if (type == 2) {
            func = cint1e_nuc_sph;
        }
    }

    int mu, nu;
    int di, dj;
    int shls[4];
    int c;
    double * buf;

    mu = 0;
    for (int i = 0; i < nbas; i++) {
        nu = 0;
        for (int j = 0; j < nbas; j++) {
            shls[0] = i; di = intcgto(i, bas);
            shls[1] = j; dj = intcgto(j, bas);

            buf = malloc(sizeof(double) * di * dj);
            func(buf, shls, atm, natm, bas, nbas, env);

            c = 0;
            for (int nuj = nu; nuj < nu + dj; nuj++) {
                for (int mui = mu; mui < mu + di; mui++) {
                    R[mui * nshells + nuj] = buf[c];
                    c++;
                }
            }
            nu += dj;
        }
        mu += di;
    }

    return R;
}

double * integral2e(int natm, int nbas, int nshells, int * atm, int * bas, double * env, int coord) {
    double * R = c_arr_doub(nshells, nshells, nshells, nshells);

    int2e func;
    cgto intcgto;

    if (coord == 0) {
        intcgto = CINTcgto_cart;
        func = cint2e_cart;
    } else if (coord == 1) {
        intcgto = CINTcgto_spheric;
        func = cint2e_sph;
    }
    
    int mu, nu, sig, lam;
    int di, dj, dk, dl;
    int shls[4];
    int c;
    double * buf;

    mu = 0;
    for (int i = 0; i < nbas; i++) {
        nu = 0;
        for (int j = 0; j < nbas; j++) {
            sig = 0;
            for (int k = 0; k < nbas; k++) {
                lam = 0;
                for (int l = 0; l < nbas; l++) {
                    shls[0] = i; di = intcgto(i, bas);
                    shls[1] = j; dj = intcgto(j, bas);
                    shls[2] = k; dk = intcgto(k, bas);
                    shls[3] = l; dl = intcgto(l, bas);

                    buf = malloc(sizeof(double) * di * dj * dk * dl);
                    func(buf, shls, atm, natm, bas, nbas, env, NULL);

                    c = 0;
                    for (int laml = lam; laml < lam + dl; laml++) {
                        for (int sigk = sig; sigk < sig + dk; sigk++) {
                            for (int nuj = nu; nuj < nu + dj; nuj++) {
                                for (int mui = mu; mui < mu + di; mui++) {
                                    R[(int) (mui*pow(nshells, 3) + nuj*pow(nshells, 2) + sigk*nshells + laml)] = buf[c];
                                    c++;
                                }
                            }
                        }
                    }

                    free(buf);

                    lam += dl;
                }
                sig += dk;
            }
            nu += dj;
        }
        mu += di;
    }

    return R;
}

void integrals(int natm, int nbas, int nshells, int * atm, int * bas, double * env, double ** S, double ** H, double ** two) 
{
    *S = c_arr(nshells, nshells);
    *H = c_arr(nshells, nshells);
    *two = c_arr_doub(nshells, nshells, nshells, nshells);

    double * T = c_arr(nshells, nshells);
    double * V = c_arr(nshells, nshells);

    int mu, nu, sig, lam;
    int di, dj, dk, dl;
    int shls[4];

    int c;
    
    double * buf;

    mu = 0;
    for (int i = 0; i < nbas; i++) {
        nu = 0;
        for (int j = 0; j < nbas; j++) {
            sig = 0;

            shls[0] = i; di = CINTcgto_cart(i, bas);
            shls[1] = j; dj = CINTcgto_cart(j, bas);

            buf = malloc(sizeof(double) * di * dj);
            cint1e_ovlp_cart(buf, shls, atm, natm, bas, nbas, env);

            c = 0;
            for (int nuj = nu; nuj < nu + dj; nuj++) {
                for (int mui = mu; mui < mu + di; mui++) {
                    (*S)[mui * nshells + nuj] = buf[c];
                    c++;
                }
            }

            c = 0;
            cint1e_kin_cart(buf, shls, atm, natm, bas, nbas, env);
            for (int nuj = nu; nuj < nu + dj; nuj++) {
                for (int mui = mu; mui < mu + di; mui++) {
                    T[mui * nshells + nuj] = buf[c];
                    c++;
                }
            }

            c = 0;
            cint1e_nuc_cart(buf, shls, atm, natm, bas, nbas, env);
            for (int nuj = nu; nuj < nu + dj; nuj++) {
                for (int mui = mu; mui < mu + di; mui++) {
                    V[mui * nshells + nuj] = buf[c];
                    c++;
                }
            }

            for (int nuj = nu; nuj < nu + dj; nuj++) {
                for (int mui = mu; mui < mu + di; mui++) {
                    (*H)[mui*nshells + nuj] = T[mui*nshells + nuj] + V[mui*nshells + nuj];
                }
            }

            free(buf);

            for (int k = 0; k < nbas; k++) {
                lam = 0;
                for (int l = 0; l < nbas; l++) {
                    shls[2] = k; dk = CINTcgto_cart(k, bas);
                    shls[3] = l; dl = CINTcgto_cart(l, bas);

                    buf = malloc(sizeof(double) * di * dj * dk * dl);
                    cint2e_cart(buf, shls, atm, natm, bas, nbas, env, NULL);

                    c = 0;
                    for (int laml = lam; laml < lam + dl; laml++) {
                        for (int sigk = sig; sigk < sig + dk; sigk++) {
                            for (int nuj = nu; nuj < nu + dj; nuj++) {
                                for (int mui = mu; mui < mu + di; mui++) {
                                    (*two)[(int) (mui*pow(nshells, 3) + nuj*pow(nshells, 2) + sigk*nshells + laml)] = buf[c];
                                    c++;
                                }
                            }
                        }
                    }

                    free(buf);

                    lam += dl;
                }
                sig += dk;
            }
            nu += dj;
        }
        mu += di;
    }
    free(T);
    free(V);
}

void find_X(int nshells, double * S, double ** X, double ** Xdag) 
{
    double * eig; // = c_arr(1, nshells);
    double * U = c_arr(nshells, nshells);
    dcopy(nshells, U, S);

    char jobz = 'V', uplo = 'U';
    integer lda = nshells, n = nshells, info = 8, lwork = 3*n;
    double w[n], work[3*n];

    dsyev_(&jobz, &uplo, &n, U, &lda, w, work, &lwork, &info);

    if(info > 0) {
        printf("The algorithm failed to compute eigenvalues.\n");
        exit(1);
    }

    double * Ut = c_arr(nshells, nshells);
    transpose(nshells, U, Ut);
    dcopy(nshells, U, Ut);

    eig = w;

    double * diag_eig = c_arr(nshells, nshells);
    for (int i = 0; i < nshells; i++) {
        for (int j = 0; j < nshells; j++) {
            if (i == j) {
                diag_eig[i*nshells + j] = pow(eig[i], -0.5);
            }
        }
    }

    *X = c_arr(nshells, nshells);
    matmult(nshells, U, diag_eig, *X);

    *Xdag = c_arr(nshells, nshells);
    transpose(nshells, *X, *Xdag);

    free(U);
    free(Ut);
    // free(eig);
    free(diag_eig);
}

void calc_F(int n, double * P, double * two, double * H, double ** F) 
{
    double * G = c_arr(n, n);
    *F = c_arr(n, n);

    for (int mu = 0; mu < n; mu++) {
        for (int nu = 0; nu < n; nu++) {
            for (int la = 0; la < n; la++) {
                for (int sig = 0; sig < n; sig++) {
                    G[mu * n + nu] += P[la * n + sig] 
                        * (two[(int) (mu*pow(n, 3) + nu*pow(n, 2) + sig*n + la)] 
                        - 0.5 * two[(int) (mu*pow(n, 3) + la*pow(n, 2) + sig*n + nu)]);
                }
            }

            (*F)[mu * n + nu] += G[mu * n + nu] + H[mu * n + nu];
        }
    }
    free(G);
}

void calc_Fprime(int n, double * F, double * X, double * Xdag, double ** Fprime) 
{
    double * inter = c_arr(n, n);
    matmult(n, Xdag, F, inter);

    *Fprime = c_arr(n, n);
    matmult(n, inter, X, *Fprime);

    free(inter);
}

void diag_F(int nshells, double * Fprime, double * X, double ** C, double ** epsilon)
{
    double * U = c_arr(nshells, nshells);
    dcopy(nshells, U, Fprime);

    char jobz = 'V', uplo = 'U';
    integer lda = nshells, n = nshells, info = 8, lwork = 3*n;
    double w[n], work[3*n];

    dsyev_(&jobz, &uplo, &n, U, &lda, w, work, &lwork, &info);

    if(info > 0) {
        printf("The algorithm failed to compute eigenvalues.\n");
        exit(1);
    }

    double * Ut = c_arr(nshells, nshells);
    transpose(nshells, U, Ut);
    dcopy(nshells, U, Ut);
    free(Ut);

    double * Udag = c_arr(nshells, nshells);
    transpose(nshells, U, Udag);
    
    double * inter = c_arr(nshells, nshells);
    matmult(nshells, Udag, Fprime, inter);
    
    double * f = c_arr(nshells, nshells);
    matmult(nshells, inter, U, f);

    double * Cprime = c_arr(nshells, nshells);
    dcopy(nshells, Cprime, U);

    *epsilon = c_arr(nshells, nshells);
    dcopy(nshells, *epsilon, f);

    *C = c_arr(nshells, nshells);
    matmult(nshells, X, Cprime, *C);
    
    free(U);
    free(Udag);
    free(inter);
    free(f);
    free(Cprime);
}

double * calc_P(int n, int nelec, double * C) 
{
    double * P = c_arr(n, n);
    for (int mu = 0; mu < n; mu++) {
        for (int nu = 0; nu < n; nu++) {
            for (int i = 0; i < (nelec / 2); i++) {
                P[mu * n + nu] += 2.0 * C[mu * n + i] * C[nu * n + i];
            }
        }
    }

    return P;
}

double f_delta(int n, double * P, double * P_old) 
{
    double delta = 0;
    for (int mu = 0; mu < n; mu++) {
        for (int nu = 0; nu < n; nu++) {
            delta += pow(P[mu * n + nu] - P_old[mu * n + nu], 2.0);
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

double * RHF(int natm, int nbas, int nelec, int nshells, int * atm, int * bas, double * env, int imax, double conv)
{   
    double * S, * H;
    double * two;
    integrals(natm, nbas, nshells, atm, bas, env, &S, &H, &two);

    double * X, * Xdag;
    find_X(nshells, S, &X, &Xdag);
    free(S);

    double * P = c_arr(nshells, nshells);
    for (int i = 0; i < nshells; i++) {
        for (int j = 0; j < nshells; j++) {
            if (i == j) {
                P[i * nshells + j] = 1.0;
            }
        }
    }

    int i = 0;
    double delta = 1.0;

    double * Pold;
    double * F;
    double * Fprime;
    double * C, * epsilon;
    
    while (delta > conv && i < imax) {
        Pold = c_arr(nshells, nshells);
        dcopy(nshells, Pold, P);
        
        calc_F(nshells, Pold, two, H, &F);
        calc_Fprime(nshells, F, X, Xdag, &Fprime);
        diag_F(nshells, Fprime, X, &C, &epsilon);

        free(P);
        P = NULL;

        P = calc_P(nshells, nelec, C);
        delta = f_delta(nshells, P, Pold);
        i++;

        free(Pold);
        Pold = NULL;
        free(F);
        F = NULL;
        free(Fprime);
        Fprime = NULL;

        free(C);
        C = NULL;
        free(epsilon);
        epsilon = NULL;

        // printf("i: %d\n", i);
        // print_arr(nshells, 2, P);
    }

    free(X);
    free(Xdag);
    free(H);
    free(two);

    if (delta > conv && i == imax) {
        printf("Did not converge\n");
        for (int i = 0; i < nshells; i++) {
            for (int j = 0; j < nshells; j++) {
                P[i*nshells + j] = 0.0;
            }
        }
    }

    // printf("conv P\n");
    // print_arr(nshells, 2, P);

    return P;
}

double energy(int natm, int nbas, int nshells, int * atm, int * bas, double * env, double * P)
{
    double * S, * H;
    double * two;
    integrals(natm, nbas, nshells, atm, bas, env, &S, &H, &two);

    double * F;
    calc_F(nshells, P, two, H, &F);

    double E0 = 0.0;
    for (int mu = 0; mu < nshells; mu++) {
        for (int nu = 0; nu < nshells; nu++) {
            E0 += 0.5 * P[mu * nshells + nu] * (H[mu * nshells + nu] + F[mu * nshells + nu]);
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

    return Enuc + E0;
}

// void energy(double * E, int natm, int nbas, int nshells, int * atm, int * bas, double * env, double * P)
// {
//     double * S, * H;
//     double * two;
//     integrals(natm, nbas, nshells, atm, bas, env, &S, &H, &two);

//     double * F;
//     calc_F(nshells, P, two, H, &F);

//     double E0 = 0.0;
//     for (int mu = 0; mu < nshells; mu++) {
//         for (int nu = 0; nu < nshells; nu++) {
//             E0 += 0.5 * P[mu * nshells + nu] * (H[mu * nshells + nu] + F[mu * nshells + nu]);
//         }
//     }

//     double Enuc = 0.0;
//     for (int i = 0; i < natm; i++) {
//         for (int j = 0; j < natm; j++) {
//             if (i > j) {
//                 Enuc += ((double) (atm[i*6 + 0] * atm[j*6 + 0]))/(norm(atm, env, i, j));
//             }
//         }
//     }

//     *E = Enuc + E0;
// }

// void energy_diff(double * E, double * dE, int natm, int nbas, int nshells, int * atm, int * bas, double * env, double * denv, double * P) {
// 	__enzyme_autodiff((void *) energy,
//         enzyme_dup, E, dE,	
//         enzyme_const, natm,
//         enzyme_const, nbas,
//         enzyme_const, nshells,
//         enzyme_const, atm,
//         enzyme_const, bas, 
//         enzyme_dup, env, denv,
//         enzyme_const, P);
// }

double * grad(int natm, int nbas, int nshells, int * atm, int * bas, double * env, double * P) 
{
    double * denv = malloc(sizeof(double) * 100);
    memset(denv, 0, sizeof(double) * 100);

    double dE = __enzyme_autodiff((void *) energy,
        enzyme_const, natm,
        enzyme_const, nbas,
        enzyme_const, nshells,
        enzyme_const, atm,
        enzyme_const, bas, 
        enzyme_dup, env, denv,
        enzyme_const, P);
    
    return denv;
}

int main()
{
    // int natm = 2;
    // int nbas = 2;
    
    // int nelec = 2;
    // int nshells = 2;

    // int * atm = malloc(sizeof(int) * natm * ATM_SLOTS);
    // int * bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
    // double * env = malloc(sizeof(double) * 100);

    // FILE * file = fopen("/u/jpmedina/libcint/molecules/h2/sto3g.txt", "r");
    // read_arrays(file, natm, nbas, &atm, &bas, &env);

    // double * ovlp = integral1e(natm, nbas, nshells, atm, bas, env, 0, 0);
    // print_arr(nshells, 2, ovlp);

    // double * two = integral2e(natm, nbas, nshells, atm, bas, env, 0);
    // print_arr(nshells, 4, two);

    // // RHF
    // int imax = 20;
    // double conv = 0.000001;

    // double * P = RHF(natm, nbas, nelec, nshells, atm, bas, env, imax, conv);

    // printf("P:\n");
    // print_arr(nshells, 2, P);

    // double E = energy(natm, nbas, nshells, atm, bas, env, P);
    // printf("E:\n%lf", E);

    // double * denv = malloc(sizeof(double) * 100);

    // double dE = __enzyme_autodiff((void *) energy,
    //     enzyme_const, natm,
    //     enzyme_const, nbas,
    //     enzyme_const, nshells,
    //     enzyme_const, atm,
    //     enzyme_const, bas, 
    //     enzyme_dup, env, denv,
    //     enzyme_const, P);

    // printf("denv:\n");
    // for (int k = 28; k < 34; k++) {
    //     printf("%f ", denv[k]);
    // }
    // printf("\n");

    // memset(denv, 0, sizeof(double) * 100);
    // printf("grad:\n");
    // denv = grad(natm, nbas, nshells, atm, bas, env, P);
    // for (int k = 28; k < 34; k++) {
    //     printf("%f ", denv[k]);
    // }
    // printf("\n");

    // double E;
    // energy(&E, natm, nbas, nshells, atm, bas, env, P);
    // printf("E: %lf\n", E);

    // double E;
    // double dE = 1.0;
    // double * denv = malloc(sizeof(double) * 10000);

    // energy_diff(&E, &dE, natm, nbas, nshells, atm, bas, env, denv, P);
    // printf("E: %lf\n", E);
    // printf("denv:\n");
    // for (int k = 28; k < 34; k++) {
    //     printf("%f ", denv[k]);
    // }
    // printf("\n");

    // printf("finite diff:\n");
    
    // double grad;
    // double E1, E2;
    // double h = 0.000001;
    
    // for (int k = 28; k < 34; k++) {
    //     env[k] += h;
    //     E1 = energy(natm, nbas, nshells, atm, bas, env, P);
    //     env[k] -= 2.0*h;
    //     E2 = energy(natm, nbas, nshells, atm, bas, env, P);
    //     env[k] += h;

    //     grad = (E1 - E2)/(2.0*h);

    //     printf("%lf ", grad);
    // }
    // printf("\n");

    return 0;
}
