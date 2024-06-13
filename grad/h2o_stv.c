#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cint.h"

int cint1e_ovlp_sph(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_kin_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

int cint1e_nuc_cart(double *buf, int *shls,
                    int *atm, int natm, int *bas, int nbas, double *env);

typedef int (* Int_cart)(double *buf, int *shls, 
					int *atm, int natm, int *bas, int nbas, double *env);

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

	// for (int i = 0; i < count; i++) {
	// 	printf("%d ", (*atm)[i]);
	// }
	// printf("\n");

	count = 0;
	while (count < nbas * BAS_SLOTS) {
        fscanf(file, "%d", &number);
		(*bas)[count] = number;
		count++;
	}

	// for (int i = 0; i < count; i++) {
	// 	printf("%d ", (*bas)[i]);
	// }
	// printf("\n");

	count = 0;
	while (fscanf(file, "%lf", &num) != EOF) {
		(*env)[count] = num;
		count++;
	}
	fclose(file);

	// for (int i = 0; i < count; i++) {
	// 	printf("%lf ", (*env)[i]);
	// }
	// printf("\n");
}

int main(int argc, char ** argv)
{	
	int natm = 3;
	int nbas = 7; // 12 for def2svp

	int * atm = malloc(sizeof(int) * natm * ATM_SLOTS);
	int * bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
	double * env = malloc(sizeof(double) * 10000);

	FILE * file = fopen("/u/jpmedina/libcint/molecules/h2o/sto3g.txt", "r");
 	read_arrays(file, natm, nbas, &atm, &bas, &env);
	
	int di, dj;
	int shls[4];

	double * buf;
	
	printf("buf:\n");

	Int_cart func = cint1e_ovlp_sph;

	for (int i = 0; i < nbas; i++) {
		for (int j = 0; j < nbas; j++) {
			shls[0] = i; di = CINTcgto_spheric(i, bas);
			shls[1] = j; dj = CINTcgto_spheric(j, bas);

			// printf("i: %d j: %d di: %d dj %d\n", i, j, di, dj);

			buf = malloc(sizeof(double) * di * dj);

			func(buf, shls, atm, natm, bas, nbas, env);

			for (int k = 0; k < di * dj; k++) {
				printf("%.3f ", buf[k]);
			}
			free(buf);
			// printf("\n");
		}
		printf("\n");
	}
	printf("\n");

	free(atm);
	free(bas);
	free(env);
}
