
#include "matUtils.h"
//#include<Accelerate/Accelerate.h>
using namespace std;

void inverse(double** A, int N)
{
    
	double *Atemp = (double*)malloc(sizeof(double)*N*N);
	int i,j;
	for(i=0; i<N; i++) {
		for(j=0; j<N; j++) {
			Atemp[j+i*N] = A[i][j];
		}
	}
	
	long int *IPIV = new long int[N+1];
	long int LWORK = N*N;
	double *WORK = new double[LWORK];
	long int INFO;
    long int Nl=N;
	dgetrf_(&Nl,&Nl,Atemp,&Nl,IPIV,&INFO);
	if(INFO != 0) {
       // std::cout<<"Error in LU-Decomposition of covariance matrix\n";
	}
    
	dgetri_(&Nl,Atemp,&Nl,IPIV,WORK,&LWORK,&INFO);
	if(INFO != 0) {
		printf("Error inverting covariance matrix.\n");
	}
	
	for(i=0; i<N; i++) {
		for(j=0; j<N; j++) {
			A[i][j] = Atemp[j+i*N];
		}
	}
	
	free(Atemp);
	delete IPIV;
	delete WORK;
}







void matrixMult(double **v1, int d11, int d12, double **v2, int d21, int d22, double **result) {
	//Rprintf("mat-mult\n");
	if(d12 != d21) {
		printf("Wrong dimensions for matrix multiplication!\n");
	}
	int i,j,k;
	
	for(i=0; i<d11; i++) {
		for(j=0; j<d22; j++) {
			result[i][j] = 0;
			for(k=0; k<d12; k++) {
				result[i][j] = result[i][j] + v1[i][k]*v2[k][j];
			}
		}
	}
}


double matrixDet(double **m, int dim) {
	    int myNCol = dim;

	    double  *myAP = new double[myNCol*(myNCol + 1)/2],
                *myW = new double[myNCol],
                *myZ = new double[myNCol*myNCol],
                *myWork = new double[myNCol * 3] ;
        long int myInfo,
        myN = (long int)(myNCol),
        myldz = (long int)(myNCol) ;

        for (register int i = 0 ; i < myN ; i++)
                for (register int j = i ; j < myldz ; j++)
                        myAP[i+(j+1)*j/2]  = m[i][j] ;

        dspev_("V", "U", &myN, myAP, myW, myZ, &myldz, myWork, &myInfo) ;

        if (myInfo != 0)
                printf("Non inversible matrix") ;
        double theDet;
        //double &theDet = theDet1;
        theDet = 1.0L ;
        for (register int i = 0 ; i < myNCol ; i++)
        {       theDet *= myW[i] ;
        }
        
        delete myAP;
        delete myW;
        delete myZ;
        delete myWork;
        
        return theDet;
}

