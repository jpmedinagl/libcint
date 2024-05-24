#include <stdio.h>
#include <stdlib.h>

#include "f2c.h"
#include "blaswrap.h"

extern int dsyev_(char *jobz, char *uplo, integer *n, doublereal *a, 
    integer *lda, doublereal *w, doublereal *work, integer *lwork, 
	integer *info);

int main() {
    char jobz='V',uplo='U';
    integer lda=3,n=3,info=8,lwork=9;
    // lapack_int lda=3,n=3,info=8;

    int i;
    double w[3], work[9];

    double a[9] = {
    3,2,4,
    2,0,2,
    4,2,3
    };


    //info=LAPACKE_dsyev(LAPACK_ROW_MAJOR,jobz,uplo,  n  ,a,  lda  , w);    
    dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);
    if( info > 0 ) {
            printf( "The algorithm failed to compute eigenvalues.\n" );
            exit( 1 );
    }
    for(i=0;i<3;i++)
    {
        printf("%f\n",w[i]);
        } 
    for(i=0;i<9;i++)
    {
            printf("%f\n",a[i]);
    }

    exit( 0 );
} 