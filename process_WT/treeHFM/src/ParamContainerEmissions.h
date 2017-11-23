
#ifndef PARAMCONTAINEREMISSIONS_HEADER
#define PARAMCONTAINEREMISSIONS_HEADER

#include "matUtils.h"
#include "MemoryAllocation.h"
#include "DebugConstants.h"
#include "EmissionFunction.h" // new added
#include <cmath>

#include <list>
using namespace std;

#define MULTIVARIATEGAUSSIAN 1
#define BERNOULLI 2
#define JOINTLYINDEPENDENT 3
#define POISSON 4
#define MULTINOMIAL 5
#define NEGATIVEBINOMIAL 6
#define POISSONLOGNORMAL 7

class ParamContainerEmissions {
  protected:
	int whichone; //!<@brief The probability density function used. 1=Gaussian. 2=Bernoulli. 3=Poisson
	int D; //!<@brief #dimensions of observation sequence.

	// multivariate gaussian
	double **mu; //!<@brief (Dx1)-dimensional matrix of emission density mean.
	double **sigma; //!<@brief (DxD)-dimensional covariance matrix of emission densities.
	double **inverseSigma; //!<@brief (DxD)-dimensional inverse covariance matrices of emission densities.
	double determinant; //!<@brief Determinant of the covariance matrix.
	double regularize; //!<@brief Regularization parameter to keep the covariance matrix invertible.
	int nsample;
	int* start;
	int* T;
	double** gammaAux;
	double logCovPrior;
	int updateCov;
	int sharedCov;

	// independent bernoulli (model introduce by Kellis et al.)
	double p;  //!<@brief Success probability of the Bernoulli trial.

	// poisson mean
	double lambda;

	// mutlinomial distribution
	double* mp;
	int n;

	// negative binomial
	double size_nb;
	double mu_nb;
	double* sizeFactor_nb;
	double pi_nb;

	// poisson log normal
	double mu_poilog;
	double sigma_poilog;
	double* sizeFactor_poilog;

	// sample x d x #unique_obs matrix
	double** uniqueObsProb;
	int** myUniqueLens;

	// more to come...
  public:
    // one constructor for each emission function.

	int getWhichOne();

      /**
	* Constructor for the Multivariate Gaussian parameter set.
	*
	* @param mu (Dx1)-dimensional matrix of emission density mean.
	* @param sigma (DxD)-dimensional covariance matrix of emission densities.
	* @param regularize Regularization parameter to keep the covariance matrix invertible.
	*
	*/
	ParamContainerEmissions(double **mu, double **sigma, double regularize, int D, int* start, int updateCov, int sharedCov); // MultivariateGaussian
	ParamContainerEmissions(double p, int D, int* start); // Bernoulli
	ParamContainerEmissions(double lambda, int D, int* start, int whichone); // Poisson
	ParamContainerEmissions(double* p, int n, int D, int* start); // Multinomial
	ParamContainerEmissions(double mu_nb, double size_nb, double* sizeFactor_nb, double pi_nb, int D, int* start); // NegativeBinomial
	ParamContainerEmissions(double mu_poilog, double sigma_poilog, double* sizeFactor_poilog, int D, int* start); // PoissonLogNormal

	ParamContainerEmissions(int D);
	double** getGaussianMU();
	double** getGaussianSIGMA();
	double** getGaussianINVSIGMA();
	double getGaussianDET();
	double getGaussianREG();
	int getSharedCov();
	double** getGammaAux();
	int* getT();
	int getUpdateCov();
	void setGaussianMU(double **mu);
	void setGaussianMUelement(double val, int d);
	void setGaussianSIGMA(double **sigma);
	void setGaussianSIGMAelement(double val, int d1, int d2);
	void setGaussianINVSIGMAelement(double val, int d1, int d2);
	void setGaussianDET(double val);
	void setDataVars(int nsample, int* T);
	void setDataVars(double** wrapper_gamma);
	void setLogCovPrior(double prior);
	double getLogCovPrior();

	void setGammaAux(double val, int n, int t);
	int getNsample();


      /**
	* Constructor for the independent Bernoulli parameter set.
	*
	* @param p D-dimensional vector of emission density mean.
	*
	*/

	double getBernoulliP();
	void setBernoulliPelement(double val, int d);
	void setBernoulliP(double p);

	double getPoissonLambda();
	void setPoissonLambdaelement(double val, int d);
	void setPoissonLambda(double lambda);

	double* getMultinomialP();
	void setMultinomialPelement(double val, int d);
	void setMultinomialP(double* p);
	int getN();

	int getD();
	int* getStart();
	virtual ~ParamContainerEmissions();

	double getMuNB();
	double getSizeNB();
	void setMuNB(double mu_nb);
	void setSizeNB(double size_nb);
	void setPiNB(double pi_nb);
	double* getSizeFactorNB();
	double getPiNB();

	double getMuPoiLog();
	double getSigmaPoiLog();
	void setMuPoiLog(double mu_poilog);
	void setSigmaPoiLog(double sigma_poilog);
	double* getSizeFactorPoiLog();


	void initUniqueObsProb(double*** observations, int nsample, int* T);
	double** getUniqueObsProb();
	int** getUniqueLens();
};
 ParamContainerEmissions** allocateParamContainerVector(int d);
#endif
