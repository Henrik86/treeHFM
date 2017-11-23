#ifndef EMISSIONFUNCTION_HEADER
#define EMISSIONFUNCTION_HEADER
 
#include "ParamContainerEmissions.h"
//#include <iostream>
#include <cmath>
#include "math.h"
#include <list>
using namespace std;

class ParamContainerEmissions;

class EmissionFunction {
	protected:
		ParamContainerEmissions *emissionParams;
		std::list<EmissionFunction*> ContainerList;
	public:

		double** updateNumeratorSIGMA;
		double** updateDenominatorSIGMA;

		virtual double calcEmissionProbability(double *obs, int isna, int currN) = 0;
		//virtual double Prior(SEXP hyperparameters) = 0;
		
		virtual void updateAuxiliariesCoupled(double*** observations, double** gamma, double* Pk, int* T, int n, int i, int statecouple, int** isNaN) = 0;
		virtual void updateAuxiliaries(double*** observations, double** gamma, double* Pk, int* T, int n, int , int** isNaNi) = 0;
		virtual void update(double ***observations, double* Pk, int** isNaN,  int currN) = 0;
		virtual ParamContainerEmissions* getParameter();
		virtual void updateCoupledRevop(double ***observations, double* Pk, int statecouple, int* state2flag, int* revop, double** revGammaAux, int** isNaN,  int currN) = 0;
		virtual void updateAuxiliariesCoupledRevop(double*** observations, double** gamma, double* Pk, int* T, int n, int i, int statecouple, int* state2flag, int* revop, int** isNaN) = 0;
		virtual void computeShared(EmissionFunction** myEmissions, int nStates);
		virtual void resetShared();
		std::list<EmissionFunction*> getEmissionFunctionList();
		//  virtual void updateParameters() = 0;
		virtual ~EmissionFunction() {}//printf("delete->EmissionFunction;\n");
};

EmissionFunction** allocateEmissionFunctionVector(int d);

#endif
