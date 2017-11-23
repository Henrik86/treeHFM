
#include "ParamContainerEmissions.h"
#include <iostream>
#define ARRAYSIZE(a) (sizeof(a) / sizeof(a[0]))

int ParamContainerEmissions::getWhichOne() {
	return this->whichone;
}

ParamContainerEmissions** allocateParamContainerVector(int d) {
	ParamContainerEmissions **vector = (ParamContainerEmissions**)malloc(sizeof(ParamContainerEmissions*)*d);
	if(vector == NULL) {
		printf("Not enough memory!\n");
	}
	return vector;
}

ParamContainerEmissions::ParamContainerEmissions(int d) {
	this->D=d;
	this->whichone = JOINTLYINDEPENDENT;

}

//ParamContainerEmissions::~ParamContainerEmissions() {}

/**
 * 01) ParamContainerEmissions Constructor for MULTIVARIATEGAUSSIAN
 *
 */
ParamContainerEmissions::ParamContainerEmissions(double **mu, double **sigma, double regularize, int D, int* start, int updateCov, int sharedCov) {
    
    //printf("create container");
    //cout<<"create container";
	this->logCovPrior = 0;
	this->updateCov = updateCov;
	this->sharedCov = sharedCov;
	this->whichone = MULTIVARIATEGAUSSIAN;
	this->mu = mu;
	this->sigma = sigma;
	this->regularize = regularize;
	this->D = D;
	this->start=start;
    
	this->inverseSigma = allocateNumericMatrix(D, D);
	int i,j;
	for(i=0; i<D; i++) {
		for(j=0; j<D; j++) {
			this->inverseSigma[i][j] = this->sigma[i][j];
            //cout<< this->inverseSigma[i][j]<<" ";

		}
        //cout<<endl;
	}
    //cout<<endl;
	inverse(this->inverseSigma, D);
	//printf("whichone: %d dimension %d , arraystartsize %d \n", this->whichone, D, ARRAYSIZE(start));
	//LapackInvAndDet(this->inverseSigma, D);

	this->determinant = matrixDet(sigma, D);
	//cout<<"determinant "<< this->determinant<<endl;
	int mem = sizeof(double*)*D + sizeof(double)*1 + 2*sizeof(double*)*D + 2*sizeof(double)*D;
	if(DEBUG_MEMORY) {
		printf("new->ParamContainerEmissions:MULTIVARIATEGAUSSIAN ; (%d bytes) ", mem);
	}
    //cout<<"new->ParamContainerEmissions:MULTIVARIATEGAUSSIAN ; (%d bytes) "<< mem<<endl;;
}
/**
 * 02) ParamContainerEmissions Constructor for BERNOULLI
 *
 */

ParamContainerEmissions::ParamContainerEmissions(double p, int D, int* start) {

	this->p = p;
	this->whichone = BERNOULLI;
	this->D = D;
	this->start=start;
	int mem = 2*sizeof(double)*D;

	if(DEBUG_MEMORY) {
			printf("new->ParamContainerEmissions:BERNOULLI; (%d bytes) ", mem);
	}
}

ParamContainerEmissions::ParamContainerEmissions(double lambda, int D, int* start, int whichone) {

	if(whichone != POISSON) {
		printf("Must be Poisson emission here!\n");
	}
	this->lambda = lambda;
	this->whichone = POISSON;
	this->D = D;
	this->start=start;
	int mem = 2*sizeof(double)*D;

	if(DEBUG_MEMORY) {
			printf("new->ParamContainerEmissions:POISSON; (%d bytes) ", mem);
	}
}

ParamContainerEmissions::ParamContainerEmissions(double* p, int n, int D, int* start) {

	this->mp = p;
	this->n = n;
	this->whichone = MULTINOMIAL;
	this->D = D;
	this->start=start;
	int mem = 2*sizeof(double)*D;

	if(DEBUG_MEMORY) {
			printf("new->ParamContainerEmissions:MULTINOMIAL; (%d bytes) ", mem);
	}
}

ParamContainerEmissions::ParamContainerEmissions(double mu_nb, double size_nb, double* sizeFactor_nb, double pi_nb, int D, int* start) {

	this->mu_nb = mu_nb;
	this->size_nb = size_nb;
	this->sizeFactor_nb = sizeFactor_nb;
	this->pi_nb = pi_nb;
	this->whichone = NEGATIVEBINOMIAL;
	this->D = D;
	this->start=start;
	//Rprintf("mu=%f, size=%f\n", this->mu_nb, this->size_nb);
	int mem = 2*sizeof(double)*D;

	if(DEBUG_MEMORY) {
			printf("new->ParamContainerEmissions:NEGATIVEBINOMIAL; (%d bytes) ", mem);
	}
}

ParamContainerEmissions::ParamContainerEmissions(double mu_poilog, double sigma_poilog, double* sizeFactor_poilog, int D, int* start) {

	this->mu_poilog = mu_poilog;
	this->sigma_poilog = sigma_poilog;
	this->sizeFactor_poilog = sizeFactor_poilog;
	this->whichone = POISSONLOGNORMAL;
	this->D = D;
	this->start=start;
	//Rprintf("mu=%f, size=%f\n", this->mu_nb, this->size_nb);
	int mem = 2*sizeof(double)*D;

	if(DEBUG_MEMORY) {
			printf("new->ParamContainerEmissions:POISSONLOGNORMAL; (%d bytes) ", mem);
	}
}

/**
 * ParamContainerEmissions DESTRUCTOR
 */
ParamContainerEmissions::~ParamContainerEmissions() {
	if(this->whichone == MULTIVARIATEGAUSSIAN) {
		int i;
		for(i=0; i<this->D; i++) {
			free(this->mu[i]);
			free(this->sigma[i]);
			free(this->inverseSigma[i]);
		}
		free(this->sigma);
        free(this->mu);
		free(this->inverseSigma);

		int mem = sizeof(double*)*this->D + sizeof(double)*1 + 2*sizeof(double*)*this->D + 2*sizeof(double)*this->D;
	
		if(DEBUG_MEMORY) {
			printf("delete->ParamContainerEmissions:MULTIVARIATEGAUSSIAN; (%d bytes) \n", mem);
		}
        //cout<<"free memory"<<endl;
		
	}
	if(this->whichone == BERNOULLI){
			int i;

			int mem = sizeof(double*)*this->D;
			if(DEBUG_MEMORY) {
				printf("delete->ParamContainerEmissions:BERNOULLI; (%d bytes) \n", mem);
			}
	}
	if(this->whichone == POISSON){
			int i;

			int mem = sizeof(double*)*this->D;
			if(DEBUG_MEMORY) {
				printf("delete->ParamContainerEmissions:BERNOULLI; (%d bytes) \n", mem);
			}
	}
    //cout<<"PME deleted"<<endl;
}


/*
 * GENERAL GETTER and SETTER
 */
void ParamContainerEmissions::setDataVars(int nsample, int* T) {
	int n,t;
	
	this->nsample = nsample;
	this->T = T;
	this->gammaAux = (double**)malloc(sizeof(double*)*nsample);
	for(n=0; n<nsample; n++) {
		this->gammaAux[n] = (double*)malloc(sizeof(double)*T[n]);
		for(t=0; t<T[n]; t++) {
			this->gammaAux[n][t] = 0;
		}
	}
}


void ParamContainerEmissions::initUniqueObsProb(double*** observations, int nsample, int* T) {
	int t, d, n, i;
	double myMax;

	this->myUniqueLens = (int**)malloc(sizeof(int*)*nsample);
	this->uniqueObsProb = (double**)malloc(sizeof(double*)*nsample);

	for(n=0; n<nsample; n++) {
		this->myUniqueLens[n] = (int*)malloc(sizeof(int)*this->D);
		for(d=0; d<this->D; d++) {
			myMax=0;
			for(t=0; t<T[n]; t++) {
				if(myMax < observations[n][t][this->start[d]]) {
					myMax = observations[n][t][this->start[d]];
				}
			}
			myMax = myMax+1;
			this->myUniqueLens[n][d] = (int)myMax;
		//	Rprintf("n=%d, d=%d, #unique=%d\n", n, this->start[d], this->myUniqueLens[n][d]);
			this->uniqueObsProb[n] = (double*)malloc(sizeof(double)*myMax);
			for(i=0; i<this->myUniqueLens[n][d]; i++) {
				this->uniqueObsProb[n][i] = -1;
			}
			for(t=0; t<T[n]; t++) {
				int currIndex = (int)(observations[n][t][this->start[d]]);
				this->uniqueObsProb[n][currIndex] = 1;
			}
		}
	}



}

double** ParamContainerEmissions::getUniqueObsProb() {
	return this->uniqueObsProb;
}

int** ParamContainerEmissions::getUniqueLens() {
	return this->myUniqueLens;
}


/*
 * GENERAL GETTER and SETTER
 */
void ParamContainerEmissions::setDataVars(double** wrapper_gamma) {
	this->gammaAux = wrapper_gamma;
}


int* ParamContainerEmissions::getT() {
	return this->T;
}


int ParamContainerEmissions::getNsample() {
	return this->nsample;
}


int ParamContainerEmissions::getD()  {
	return this->D;
}


int* ParamContainerEmissions::getStart()  {
	return this->start;
}

/*
 * MULTIVARIATEGAUSSIAN GETTER and SETTER
 */

void ParamContainerEmissions::setGammaAux(double val, int n, int t) {
	this->gammaAux[n][t] = val;
}
double** ParamContainerEmissions::getGammaAux() {
	return this->gammaAux;
}

double** ParamContainerEmissions::getGaussianMU() {
	return this->mu;
}
double** ParamContainerEmissions::getGaussianSIGMA() {
	return this->sigma;
}
double** ParamContainerEmissions::getGaussianINVSIGMA() {
	return this->inverseSigma;
}
double ParamContainerEmissions::getGaussianDET() {
	return this->determinant;
}
double ParamContainerEmissions::getGaussianREG()  {
	return this->regularize;
}

int ParamContainerEmissions::getUpdateCov() {
	return this->updateCov;
}

void ParamContainerEmissions::setLogCovPrior(double prior) {
	this->logCovPrior = prior;
}

double ParamContainerEmissions::getLogCovPrior() {
	return this->logCovPrior;
}

void ParamContainerEmissions::setGaussianMU(double **mu) {
	int i;
	for(i=0; i<D; i++) {
		this->mu[i][0] = mu[i][0];
	}	
}

void ParamContainerEmissions::setGaussianMUelement(double val, int d) {
	this->mu[d][0] = val;
}

void ParamContainerEmissions::setGaussianSIGMAelement(double val, int d1, int d2) {
	this->sigma[d1][d2] = val;
}

void ParamContainerEmissions::setGaussianINVSIGMAelement(double val, int d1, int d2) {
	this->inverseSigma[d1][d2] = val;
}

void ParamContainerEmissions::setGaussianDET(double val) {
	this->determinant = val;
}

void ParamContainerEmissions::setGaussianSIGMA(double **sigma) {
	int i,j;
	for(i=0; i<D; i++) {
		for(j=0; j<D; j++) {
			this->sigma[i][j] = sigma[i][j];
			this->inverseSigma[i][j] = sigma[i][j];
		}
	}
	inverse(inverseSigma, this->D);
	this->determinant = matrixDet(sigma, this->D);
}


/*
 * BERNOULLI GETTER and SETTER
 *
 */

double ParamContainerEmissions::getBernoulliP(){
	return this->p;
}

void ParamContainerEmissions::setBernoulliPelement(double val, int d){
	this->p = val;

}

void ParamContainerEmissions::setBernoulliP(double p){
	this->p = p;

}


double  ParamContainerEmissions::getPoissonLambda() {
	return this->lambda;
}

void  ParamContainerEmissions::setPoissonLambdaelement(double val, int d) {
	this->lambda = val;
}
void  ParamContainerEmissions::setPoissonLambda(double lambda) {
	this->lambda = lambda;
}


double* ParamContainerEmissions::getMultinomialP(){
	return this->mp;
}

int ParamContainerEmissions::getN(){
	return this->n;
}


void ParamContainerEmissions::setMultinomialPelement(double val, int d) {
	this->mp[d] = val;
}


void ParamContainerEmissions::setMultinomialP(double* p) {
	this->mp = p;
}



double ParamContainerEmissions::getMuNB() {
	return this->mu_nb;
}

double ParamContainerEmissions::getSizeNB() {
	return this->size_nb;
}
double ParamContainerEmissions::getPiNB() {
	return this->pi_nb;
}

void ParamContainerEmissions::setMuNB(double mu_nb) {
	this->mu_nb = mu_nb;
}
void ParamContainerEmissions::setSizeNB(double size_nb) {
	this->size_nb = size_nb;
}
void ParamContainerEmissions::setPiNB(double pi_nb) {
	this->pi_nb = pi_nb;
}

double* ParamContainerEmissions::getSizeFactorNB() {
	return this->sizeFactor_nb;
}

int ParamContainerEmissions::getSharedCov() {
	return this->sharedCov;
}


double ParamContainerEmissions::getMuPoiLog() {
	return this->mu_poilog;
}
double ParamContainerEmissions::getSigmaPoiLog() {
	return this->sigma_poilog;
}
void ParamContainerEmissions::setMuPoiLog(double mu_poilog) {
	this->mu_poilog = mu_poilog;
}
void ParamContainerEmissions::setSigmaPoiLog(double sigma_poilog) {
	this->sigma_poilog = sigma_poilog;
}
double* ParamContainerEmissions::getSizeFactorPoiLog() {
	return this->sizeFactor_poilog;
}
