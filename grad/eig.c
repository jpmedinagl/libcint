#include <stdio.h>
#include <stdlib.h>

#include "f2c.h"
#include "blaswrap.h"

extern int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a, 
    integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *info);

int main() {
    char jobz = 'V', uplo = 'U';
    integer lda = 2, n = 2, info = 8, lwork = 6;

    int i;
    double w[2], work[6];

    double a[4] = {
        1.000000, 0.222082,
        0.222082, 1.000000
    };

    dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
    if(info > 0) {
        printf("The algorithm failed to compute eigenvalues.\n");
        return 1;
    }

    for(i = 0; i < 2; i++) {
        printf("%f ", w[i]);
    }
    printf("\n");

    for(i = 0; i < 4; i++) {
        printf("%f ", a[i]);
    }
    printf("\n");

    return 0;
} 