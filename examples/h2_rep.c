#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cint.h"


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

int main(int argc, char ** argv)
{	
	int natm = 2;
	int nbas = 2;

	int * atm = malloc(sizeof(int) * natm * ATM_SLOTS);
	int * bas = malloc(sizeof(int) * nbas * BAS_SLOTS);
	double * env = malloc(sizeof(double) * 10000);

	FILE * file = fopen("/u/jpmedina/libcint/molecules/h2/sto3g.txt", "r");
 	read_arrays(file, natm, nbas, &atm, &bas, &env);
	
	int di, dj, dk, dl;
	int shls[4];

	double * buf;
	
	printf("buf:\n");

	for (int i = 0; i < nbas; i++) {
		for (int j = 0; j < nbas; j++) {
			for (int k = 0; k < nbas; k++) {
				for (int l = 0; l < nbas; l++) {
					shls[0] = i; di = CINTcgto_cart(i, bas);
					shls[1] = j; dj = CINTcgto_cart(j, bas);
					shls[2] = k; dk = CINTcgto_cart(k, bas);
					shls[3] = l; dl = CINTcgto_cart(l, bas);

					buf = malloc(sizeof(double) * di * dj * dk * dl);
					cint2e_cart(buf, shls, atm, natm, bas, nbas, env, NULL);

					for (int m = 0; m < di * dj * dk * dl; m++) {
						printf("%f ", buf[m]);
					}
					free(buf);
				}
				printf("\n");
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("\n");

	free(atm);
	free(bas);
	free(env);
}
