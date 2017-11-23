

#ifndef FGBMA_BaumWelch_h
#define FGBMA_BaumWelch_h

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include "FactorGraph.h"
#include "SumProduct.h"
#include "helpFunctions.h"
#include "MultivariateGaussian.h"

struct observation{
    std::vector<std::vector<double> > data;
};
struct observationSequence{
    std::vector<std::vector<double> > data;
};
struct observationSequenceDiscrete{
    std::vector<int> data;
};
struct allObservationSequence{
    std::vector<observationSequence> data;
};
struct allObservationSequenceDiscrete{
    std::vector<observationSequence> data;
};
struct eijDiv{
    std::vector<std::vector<std::vector<double> > > data;
};
struct allEijDiv{
    std::vector<eijDiv> data;
};
struct eijSeq{
    std::vector<std::vector<std::vector<double> > > data;
};
struct allEijSeq{
    std::vector<eijSeq> data;
};
struct ypsilon{
    std::vector<std::vector<double> > data;
};
struct parentIndices{
    std::vector<int>  data;
};
struct nodeIndices{
    std::vector<int>  data;
};
struct allYpsilon{
    std::vector<ypsilon> data;
};
struct outputC{
    std::vector<std::vector<double> > transProbSeq;
    std::vector<std::vector<double> > transProbDiv;
    std::vector<double> allLL;
    std::vector<std::vector<double> > allLLSingle;
    std::vector<std::vector<double> > mu;
    std::vector<std::vector<std::vector<double> > > sigma;
    std::vector<std::vector<double> > prior;
    
};
struct outputD{
    std::vector<std::vector<double> > transProbSeq;
    std::vector<std::vector<double> > transProbDiv;
    std::vector<std::vector<double> > prior;
    std::vector<double> allLL;
    std::vector<std::vector<double> > allLLSingle;
    std::vector<std::vector<double> > emProb;
    
};
std::vector<std::vector<std::vector<double> > > calculateEijSeq(factorGraph &fgS,std::vector<std::vector<double> >  &transProbSeq,std::vector<std::vector<double> > &observationProbabilities,int numberNodes );
std::vector<std::vector<std::vector<double> > >  calculateEijDiv(factorGraph &fgS,std::vector<std::vector<double> >  &transProbDivN,std::vector<std::vector<double> > &observationProbabilities,std::vector<std::vector<int> > stateIndicesSingle, int numberNodes, int nHStates, int nDStates);
std::vector<std::vector<double> > calculateYpsilon(factorGraph &fgS,std::vector<std::vector<double> > &observationProbabilities,int numberNodes);
std::vector<std::vector<std::vector<double> > >  calculateSigma(std::vector<observationSequence> allObservationSequences, int nHStates, std::vector<ypsilon> allYpsilons,std::vector <int> skipSequences,std::vector<std::vector<double> > muN,double covarianceFaktor, int dimData,std::vector<std::vector<std::vector<double> > >  covOld );
std::vector<double> calculatePI(int nHStates,std::vector<ypsilon> allYpsilon, std::vector <int> skipSequences);
std::vector<std::vector<double> > calculateTransProbSeq(int nHStates, std::vector <int> skipSequences,std::vector <eijSeq>  allEijSeq );
std::vector<std::vector<double> > calculateTransProbDiv(int nHStates,int nDStates, std::vector <int> skipSequences,std::vector <eijDiv>  allEijDiv);
std::vector<std::vector<double> > calculateMu(std::vector <observationSequence> allObservationSequences,int nHStates, std::vector <int> skipSequences,std::vector <ypsilon> allYpsilon, int dimData);
outputC baumWelch( std::vector <observationSequence> allObservationSequences, std::vector<std::vector<double> > startProbN, std::vector<std::vector<double> > transProbSeqN, std::vector<std::vector<double> > transProbDivN, std::vector<std::vector<double> > emProbN, std::vector<std::vector<double> > muN,  std::vector<std::vector<std::vector<double> > > SigmaN, std::vector<nodeIndices> allNodeIndices, std::vector<std::vector<int> > summationIndices, std::vector<std::vector<int> > stateIndicesSingle, std::vector<int> indicesXD, std::vector<parentIndices > allParentIndices, int type);
outputD baumWelchDiscrete( std::vector <observationSequenceDiscrete> allObservationSequences, std::vector<std::vector<double> > startProbN, std::vector<std::vector<double> > transProbSeqN, std::vector<std::vector<double> > transProbDivN, std::vector<std::vector<double> > emProbN,  std::vector<nodeIndices> allNodeIndices, std::vector<std::vector<int> > summationIndices, std::vector<std::vector<int> > stateIndicesSingle, std::vector<int> indicesXD, std::vector<parentIndices > allParentIndices, int type);
#endif
