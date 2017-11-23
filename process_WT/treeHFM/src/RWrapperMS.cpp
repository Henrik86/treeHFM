//
//  MatlabWrapperMS.cpp
//  FGBMA
//
//  Created by Henrik Failmezger on 16.01.15.
//  Copyright (c) 2015 Henrik Failmezger. All rights reserved.
//
#include "BaumWelch.h"
#include "RWrapperFunctions.h"
#include "MaxSum.h"



using namespace std;
extern "C" {
SEXP RWrapperMS(SEXP sexpstartProbN, SEXP  sexptransProbSeqN, SEXP sexptransProbDiv, SEXP sexpemProbN, SEXP sexpobservationSequences, SEXP sexpnodeIndices,SEXP sexpparents, SEXP sexptype, SEXP sexpnHStates,SEXP sexpnDStates, SEXP sexplenObs, SEXP sexpdimData, SEXP sexpsigma, SEXP sexpmu, SEXP sexpsummationIndices, SEXP sexpindicesXD, SEXP sexpnEStates){
   
    int nHStates=INTEGER(sexpnHStates)[0];
    //cout<<"nDStates"<<endl;
    int nDStates=INTEGER(sexpnDStates)[0];
    //cout<<"lenObs"<<endl;
    int lenObs=INTEGER(sexplenObs)[0];
    //cout<<"lenObs: "<<lenObs<<endl;
    int nEStates=INTEGER(sexpnEStates)[0];
    int dimData=INTEGER(sexpdimData)[0];
    //cout<<"dimData: "<<dimData<<endl;
    //Rprintf("Type\n");
    const char *typeC = CHAR(STRING_ELT(sexptype, 0));
    char type=typeC[0];
    //////////////////////// TYPE //////////////////////////////////////////////
    ////////////////////////////////////////////////////////
    
    //transProbSeqM////////////////////////////////////////////////////////
    std::vector<std::vector<double> > startProbN=  RGETMAT(sexpstartProbN, nHStates, nHStates);
    //cout<<"startProbN"<<endl;
    //printVector(startProbN);
    //cout<<"startProbN"<<endl;
    //printVector(startProbN);
    ///////////////////////////////////////////////////////////////////////
    std::vector<std::vector<double> > transProbSeqN= RGETMAT(sexptransProbSeqN, nHStates, nHStates);
    //////////transProbDiv/////////////////////////////////////////////////
    std::vector<std::vector<double> > transProbDiv= RGETMAT(sexptransProbDiv, nHStates, nDStates);
    ///////////////////////////////////////////////////////////////////////
    ///////////probO///////////////////////////////////////////////////////
    //cout<<"nEStates: "<<nEStates<<endl;
    std::vector<std::vector <double> > emProbN= RGETMAT(sexpemProbN, lenObs, nEStates);
    ///////////////////////////////////////////////////////////////////////
    
    //cout<<"observation sequences"<<endl;
    ////////////observation sequences//////////////////////////////////////
    ////////////observation sequences//////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////
    
    ////////////node indices//////////////////////////////////////
    std::vector <int> nodeIndices= RGETVECTINT(sexpnodeIndices, lenObs);
    ////////////parents//////////////////////////////////////
    //cout<<"parents: "<<endl;
    std::vector <int> parents= RGETVECTINT(sexpparents, lenObs);
    ///////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////
    ////////////summation indices//////////////////////////////////////
    //cout<<"summation indices: "<<endl;
    std::vector <std::vector <int> > summationIndices= RGETMATINT(sexpsummationIndices, nHStates, nHStates);
    summationIndices=subEl(summationIndices,1);
    printVector(summationIndices);
    ///////////////////////////////////////////////////////////////////////
    ////////////indices XD//////////////////////////////////////
    R_len_t IXD = length(sexpindicesXD);
    //INTEGER(RdimXD)[1];
    //cout<<"indices XD: "<<IXD<<endl;
    std::vector <int> indicesXD=RGETVECTINT(sexpindicesXD, IXD);
    indicesXD= subEl(indicesXD,1); //bring them to c++ indices
    //////////////////////// StateIndicesSingle //////////////////////////////////////////////

    //////////////////////// SigmaInit //////////////////////////////////////////////
    //cout<<"Sigma : "<<endl;
    std::vector<std::vector<std::vector<double> > > SigmaInit=RGETMAT(sexpsigma, dimData, nHStates, nHStates);
    //
    ///////////////////////////////////////////////////////////////////////
    //////////////////////// MuInit //////////////////////////////////////////////
    std::vector<std::vector<double> > muN;
    if(type=='c'){
        muN=RGETMAT(sexpmu, nHStates, dimData);
    }
    ///////////////////////////////////////////////////////////////////////

    int vSize=transProbSeqN.size();
    
    std::vector< std::vector<std::vector <double> > > transProbDivUD;
    transProbDivUD.push_back(transProbDiv);
    
    transProbDivUD.push_back(arrangeItems(transProbDiv,indicesXD));
    
    transitionMatrices tM;
    tM.probDiv=transProbDivUD;
    std::vector< std::vector<std::vector<double> > > startProbT;
    startProbT.push_back(startProbN);
    tM.probStart=startProbT;
    std::vector< std::vector<std::vector<double> > > transProbSeqT;
    transProbSeqT.push_back(transProbSeqN);
    tM.probSeq=transProbSeqT;
    factorGraph fg;
    //////////// CALCULATE OBSERVATION PROBABILITIES ///////////////////////////////////
    int obsLen;
    
    std::vector<std::vector<double> > probO(lenObs,std::vector<double>(nHStates));
    if(type=='d'){
        //std::cout<<"observationSequences"<<std::endl;
        std::vector <int> observationSequencesINT=RGETVECTINT(sexpobservationSequences, lenObs);
        obsLen=observationSequencesINT.size();
        //
        //cout<<"lenObs: "<<lenObs<<endl;
        //cout<<"observationSequencesINT: "<<observationSequencesINT.size()<<endl;
        //printVector(observationSequencesINT);
        //printVector(emProbN);
        for(int d=0;d< observationSequencesINT.size(); d++){
            for(int s=0; s < nHStates; s++){
                //cout<<"d: "<<d<<"s: "<<s<<endl;
                //cout<<emProbN.at(s).at(observationSequencesINT.at(d)-1)<<endl;
                probO.at(d).at(s)=emProbN.at(s).at(observationSequencesINT.at(d)-1);
            }
        }
    }else{
        SEXP Rdim=getAttrib(sexpobservationSequences,R_DimSymbol);
        int IOS=INTEGER(Rdim)[0];
        int JOS=INTEGER(Rdim)[1];
        observationSequence obsSequence;
        //cout<<"Rdim: IOS: "<<IOS<<"JOS: "<<JOS<<endl;
        std::vector<std::vector <double> >  observationSequences= RGETMAT(sexpobservationSequences, IOS, JOS);
        obsLen=JOS;
        int dimensionData=observationSequences.size();
        int* start_d = (int*)malloc(sizeof(int)*dimensionData);
        int o;
        for (o = 0; o < dimensionData; o++) {
            start_d[o] = o;
            //Rprintf("init %d \n", start_d[o]);
        }
        //
        for(int d=0;d< obsLen; d++){
            for(int s=0; s < nHStates; s++){
                ParamContainerEmissions* multGParams=new ParamContainerEmissions(vectorToArray(getColumn(muN,s)), vectorToArray(SigmaInit.at(s)), 0, dimensionData, start_d, 1, 1);
                MultivariateGaussian mG=MultivariateGaussian(multGParams);
                float prob = mG.calcEmissionProbability(vectorToArray1D(getColumn(observationSequences,d)), 0, 0);
                //cout<<"prob: "<<prob<<endl;
                probO.at(d).at(s)=prob;
            }
        }
        //
        free (start_d);
    }
    //////////////////////////////////////////////////////////////////////
    createHMT(fg,nodeIndices, parents,summationIndices,1,vSize,indicesXD,type,probO,transProbSeqN.size(),transProbDiv.at(0).size());
    std::vector<int> startNodes;
    startNodes.push_back(1);
    std::vector<std::vector<int> >maxPath=maxSumAlgorithm(fg,fg.endNodes,startNodes,tM,obsLen);
    //cout <<"Time for SumProduct (seconds): "<< ((double)(finish - start))/CLOCKS_PER_SEC<<endl;
    ///////////////////////////////////////////////////////////////////////
    SEXP Viterbi, wnames;
    PROTECT(Viterbi = NEW_LIST(3));
    PROTECT(wnames = NEW_CHARACTER(3));
    SET_STRING_ELT(wnames, 0, mkChar("maxPath"));
    SET_STRING_ELT(wnames, 1, mkChar("strainParents"));
    SET_STRING_ELT(wnames, 2, mkChar("dims"));
    SET_NAMES(Viterbi, wnames);
    UNPROTECT(1);
    SEXP maxPathOut = RPREPAREMAT(maxPath);
    SET_ELEMENT(Viterbi, 0, maxPathOut);
    ///
    ///
    //Rprintf("transProbDiv: \n");
    //printVector(fg.strainParents);
    SEXP strainParentsOut = RPREPAREMAT(fg.strainParents);
    SET_ELEMENT(Viterbi, 1, strainParentsOut);
    ///
    std::vector <int> dims;
    dims.push_back(fg.strainParents.size());
    dims.push_back(fg.strainParents.at(0).size());
    //Rprintf("transProbDiv: \n");
    //printVector(fg.strainParents);
    SEXP dimsOut = RPREPAREVECT(dims);
    SET_ELEMENT(Viterbi, 2, dimsOut);
    UNPROTECT(1);
    ///////////////////////////////////////////////////////////////////////
        //for (int i=1; i<(fgN.factorNodes.size()); i++) {
    //
    //}
    //cout <<"Time for Pass out(seconds): "<< ((double)(finishPass - startPass))/CLOCKS_PER_SEC<<endl;
    //cout<<"--------------------------------------"<<endl;
    return Viterbi;
}
}
