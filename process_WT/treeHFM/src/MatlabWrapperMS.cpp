//
//  MatlabWrapperMS.cpp
//  FGBMA
//
//  Created by Henrik Failmezger on 16.01.15.
//  Copyright (c) 2015 Henrik Failmezger. All rights reserved.
//

#include "MatlabWrapperMS.h"
#include "MaxSum.h"
#include "mex.h"
#include "MultivariateGaussian.h"

using namespace std;

mxArray * getMexArray(std::vector<double> v){
    
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    return mx;
}

mxArray * getMexArrayDeleteFE(std::vector<double> v){
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    for (int i=1;i<(v.size());i++){ //starts with one!!
        mxGetPr(mx)[i-1]=v.at(i);
    }
    return mx;
}
mxArray * getMexArray(int v){
    mxArray * mx = mxCreateDoubleMatrix(1,1, mxREAL);
    mxGetPr(mx)[0]=v;
    return mx;
}
mxArray * getMexArray(double v){
    
    mxArray * mx = mxCreateDoubleMatrix(1,1, mxREAL);
    mxGetPr(mx)[0]=v;
    return mx;
}
mxArray * getMexArray(std::vector<int> v){
    mxArray * mx = mxCreateDoubleMatrix(1,v.size(), mxREAL);
    std::copy(v.begin(), v.end(), mxGetPr(mx));
    return mx;
}
mxArray * getMexArray(std::vector<std::vector<double> >v){
    mxArray * mx = mxCreateDoubleMatrix(v.at(0).size(),v.size(), mxREAL);
    for (int i=0;i<v.size();i++){
        for (int j=0;j<v.at(i).size();j++){
            //std::copy(v.at(i).begin(), v.at(i).end(), arr);
            mxGetPr(mx)[i+(j*v.size())]=v.at(i).at(j);
        }
    }
    return mx;
}
mxArray * getMexArray(std::vector<std::vector<int> >v){
    mxArray * mx;
    if(v.size()>0){
        mx = mxCreateDoubleMatrix(v.size(),v.at(0).size(), mxREAL);
        int counter=0;
        if(v.size()>0){
            for (int i=0;i<v.at(0).size();i++){
                for (int j=0;j<v.size();j++){
                    //double arr[v.at(i).size()];
                    //std::copy(v.at(i).begin(), v.at(i).end(), arr);
                    mxGetPr(mx)[counter]=v.at(j).at(i);
                    counter++;
                }
            }
        }
    }else{
        mx = mxCreateDoubleMatrix(0,0, mxREAL);
    }
    
    return mx;
}
mxArray * createCellArray(std::vector<std::vector<double> > v, int sizeMessages){
    mxArray * mx = mxCreateCellMatrix(1,v.size());
    for (int i=1;i<v.size();i++){
        mxSetCell(mx,i-1,getMexArray(v.at(i)));
    }
    return mx;
}

mxArray * createCellArray(std::vector<std::vector<std::vector<double> > >v, int sizeMessages){
    mxArray * mx = mxCreateCellMatrix(1,v.size());
    for (int i=1;i<v.size();i++){
        mxSetCell(mx,i-1,getMexArray(v.at(i)));
    }
    return mx;
}

void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] ){
   
    clock_t startLoadData, finishLoadData;
    startLoadData=clock();
    //
    //////////////////////// TYPE //////////////////////////////////////////////
    mxChar * pData9 = mxGetChars(prhs[9]);
    char type = pData9[0];
    //startProb transProbSeqN transProbDiv emProbN nodeIndices parents summationIndices indicesXD stateIndicesSingle type lenObs mu sigma observationSequenceC observationSequenceD dimensionData
    //////////////////////// TYPE //////////////////////////////////////////////

    
    ////////////node indices//////////////////////////////////////
    mwSize nDims4 = mxGetNumberOfDimensions(prhs[4]);
    double * pData4 = mxGetPr(prhs[4]);
    std::vector <int> nodeIndices(mxGetDimensions(prhs[4])[1]);
    for(int i=0;i<mxGetDimensions(prhs[4])[1];i++){
        nodeIndices.at(i)=pData4[i];
    }
    ////////////parents//////////////////////////////////////
    mwSize nDims5 = mxGetNumberOfDimensions(prhs[5]);
    double * pData5 = mxGetPr(prhs[5]);
    std::vector <int> parents(mxGetDimensions(prhs[5])[1]);
    for(int i=0;i<mxGetDimensions(prhs[5])[1];i++){
        parents.at(i)=pData5[i];
    }
    ///////////////////////////////////////////////////////////////////////
    int lenObs=nodeIndices.size();
   
    //////////////////////////////////////////////
    double *sN;
    sN= mxGetPr (prhs[0]);
    mwSize n = mxGetNumberOfElements (prhs[0]);
    //transProbSeqM////////////////////////////////////////////////////////
    mwSize nDims0 = mxGetNumberOfDimensions(prhs[0]);
    double * pData0 = mxGetPr(prhs[0]);
    std::vector<std::vector<double> > startProbN(mxGetDimensions(prhs[0])[0], std::vector<double>(mxGetDimensions(prhs[0])[1]));
    for(int i=0;i<mxGetDimensions(prhs[0])[0];i++){
        for(int j=0;j<mxGetDimensions(prhs[0])[1];j++){
            startProbN.at(i).at(j)=pData0[i+j*mxGetDimensions(prhs[0])[0]];
        }
    }
    ///////////////////////////////////////////////////////////////////////
    mwSize nDims1 = mxGetNumberOfDimensions(prhs[1]);
    double * pData1 = mxGetPr(prhs[1]);
    std::vector<std::vector<double> > transProbSeqN(mxGetDimensions(prhs[1])[0], std::vector<double>(mxGetDimensions(prhs[1])[1]));
    for(int i=0;i<mxGetDimensions(prhs[1])[0];i++){
        for(int j=0;j<mxGetDimensions(prhs[1])[1];j++){
            transProbSeqN.at(i).at(j)=pData1[i+j*mxGetDimensions(prhs[1])[0]];
        }
    }
    int nHStates=mxGetDimensions(prhs[1])[0];
    //////////transProbDiv/////////////////////////////////////////////////
    mwSize nDims2 = mxGetNumberOfDimensions(prhs[2]);
    double * pData2 = mxGetPr(prhs[2]);
    std::vector<std::vector<double> > transProbDiv(mxGetDimensions(prhs[2])[0], std::vector<double>(mxGetDimensions(prhs[2])[1]));
    for(int i=0;i<mxGetDimensions(prhs[2])[0];i++){
        for(int j=0;j<mxGetDimensions(prhs[2])[1];j++){
            transProbDiv.at(i).at(j)=pData2[i+j*mxGetDimensions(prhs[2])[0]];
        }
    }
    ///////////////////////////////////////////////////////////////////////
    ///////////probO///////////////////////////////////////////////////////
    mwSize nDims3 = mxGetNumberOfDimensions(prhs[3]);
    double * pData3 = mxGetPr(prhs[3]);
    std::vector<std::vector <double> > emProbN(mxGetDimensions(prhs[3])[0], std::vector<double>(mxGetDimensions(prhs[3])[1]));
    for(int i=0;i<mxGetDimensions(prhs[3])[0];i++){
        for(int j=0;j<mxGetDimensions(prhs[3])[1];j++){
            emProbN.at(i).at(j)=pData3[i+j*mxGetDimensions(prhs[3])[0]];
        }
    }

    ///////////////////////////////////////////////////////////////////////
    
    ///////////probO///////////////////////////////////////////////////////
   
    std::vector<std::vector<double> > probO(lenObs,std::vector<double>(nHStates));
    if(type=='d'){
 
        mwSize nDims8 = mxGetNumberOfDimensions(prhs[13]);
        double * pData8 = mxGetPr(prhs[13]);
        std::vector <int> observationSequencesINT(lenObs);
        for(int i=0;i<mxGetDimensions(prhs[13])[1];i++){
            observationSequencesINT.at(i)=pData8[i];
        }
        

        for(int d=0;d< observationSequencesINT.size(); d++){
            for(int s=0; s < nHStates; s++){
                probO.at(d).at(s)=emProbN.at(s).at(observationSequencesINT.at(d)-1);
            }
        }
    }else{
        ////

        mwSize nDims11 = mxGetNumberOfDimensions(prhs[11]);
        double * pData11 = mxGetPr(prhs[11]);
        std::vector<std::vector<std::vector<double> > > SigmaInit(mxGetDimensions(prhs[11])[2], std::vector< std::vector<double> >(mxGetDimensions(prhs[11])[1], std::vector<double> (mxGetDimensions(prhs[11])[0])));
        //
        
        //
        for(int k=0;k<mxGetDimensions(prhs[11])[2];k++){
            for(int i=0;i<mxGetDimensions(prhs[11])[0];i++){
                for(int j=0;j<mxGetDimensions(prhs[11])[1];j++){
                    SigmaInit.at(k).at(i).at(j)=pData11[i+j*mxGetDimensions(prhs[11])[0]+k*mxGetDimensions(prhs[11])[0]* mxGetDimensions(prhs[11])[1]];
                }
            }
        }
        ///////////////////////////////////////////////////////////////////////
        //////////////////////// MuInit //////////////////////////////////////////////
        mwSize nDims12 = mxGetNumberOfDimensions(prhs[12]);
        double * pData12 = mxGetPr(prhs[12]);
        std::vector<std::vector<double> > muN(mxGetDimensions(prhs[12])[0], std::vector<double>(mxGetDimensions(prhs[12])[1]));
        //
        for(int i=0;i<mxGetDimensions(prhs[12])[0];i++){
            for(int j=0;j<mxGetDimensions(prhs[12])[1];j++){
                muN.at(i).at(j)=pData12[i+j*mxGetDimensions(prhs[12])[0]];
            }
        }
        ////////////////////////
        mwSize nDims13 = mxGetNumberOfDimensions(prhs[13]);
        double * pData13 = mxGetPr(prhs[13]);
        std::vector<std::vector <double> > observationSequence(mxGetDimensions(prhs[13])[0], std::vector<double>(mxGetDimensions(prhs[13])[1]));
        for(int i=0;i<mxGetDimensions(prhs[13])[0];i++){
            for(int j=0;j<mxGetDimensions(prhs[13])[1];j++){
                observationSequence.at(i).at(j)=pData13[i+j*mxGetDimensions(prhs[13])[0]];
            }
        }
        ////////////////////////
        int dimensionData=observationSequence.size();
        ////////////////////////
        int* start_d = (int*)malloc(sizeof(int)*dimensionData);
        int o;
        for (o = 0; o < dimensionData; o++) {
            start_d[o] = o;
        }
        //
        for(int d=0;d< lenObs; d++){
            for(int s=0; s < nHStates; s++){
                ParamContainerEmissions* multGParams=new ParamContainerEmissions(vectorToArray(getColumn(muN,s)), vectorToArray(SigmaInit.at(s)), 0, dimensionData, start_d, 1, 1);
                MultivariateGaussian* mG= new MultivariateGaussian(multGParams);
                float prob = mG->calcEmissionProbability(vectorToArray1D(getColumn(observationSequence,d)), 0, 0);
                probO.at(d).at(s)=prob;
                delete mG;
            }
        }
        //
        free (start_d);
    }
    ///////////////////////////////////////////////////////////////////////
    ////////////observation sequences//////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////
    ////////////summation indices//////////////////////////////////////
    mwSize nDims6 = mxGetNumberOfDimensions(prhs[6]);
    double * pData6 = mxGetPr(prhs[6]);
    std::vector <std::vector <int> > summationIndices(mxGetDimensions(prhs[6])[0], std::vector<int>(mxGetDimensions(prhs[6])[1]));
    for(int i=0;i<mxGetDimensions(prhs[6])[0];i++){
        for(int j=0;j<mxGetDimensions(prhs[6])[1];j++){
            summationIndices.at(i).at(j)=pData6[i+j*mxGetDimensions(prhs[6])[0]];
        }
    }
    summationIndices=subEl(summationIndices,1);
    ///////////////////////////////////////////////////////////////////////
    ////////////indices XD//////////////////////////////////////
    mwSize nDims7 = mxGetNumberOfDimensions(prhs[7]);
    double * pData7 = mxGetPr(prhs[7]);
    std::vector <int> indicesXD(mxGetDimensions(prhs[7])[1]);
    for(int i=0;i<mxGetDimensions(prhs[7])[1];i++){
        indicesXD.at(i)=pData7[i];
    }
    indicesXD= subEl(indicesXD,1); //bring them to c++ indices
    //////////////////////// StateIndicesSingle //////////////////////////////////////////////
    mwSize nDims8 = mxGetNumberOfDimensions(prhs[8]);
    double * pData8 = mxGetPr(prhs[8]);
    std::vector<std::vector <int> > stateIndicesSingle(mxGetDimensions(prhs[8])[0], std::vector<int>(mxGetDimensions(prhs[8])[1]));
    for(int i=0;i<mxGetDimensions(prhs[8])[0];i++){
        for(int j=0;j<mxGetDimensions(prhs[8])[1];j++){
            stateIndicesSingle.at(i).at(j)=pData8[i+j*mxGetDimensions(prhs[8])[0]];
        }
    }
    finishLoadData=clock();
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
    //
    
    factorGraph fg;
    createHMT(fg,nodeIndices, parents,summationIndices,1,vSize,indicesXD,type,probO,transProbSeqN.size(),transProbDiv.at(0).size());
   //////////////////////////
    
   //////////////////////////
    std::vector<int> startNodes;
    startNodes.push_back(1);
    
    int topologyCheck=1;
    if(fg.strainParents.size()>0){
        std::vector<double>  tP=tabulate(getColumn(fg.strainParents,0));
        Max maxParents;
        maxParents = max(tP);
        if(maxParents.value >=3){
            topologyCheck=0;
        }
    }
    if(topologyCheck==1){
        std::vector<std::vector<int> >maxPath=maxSumAlgorithm(fg,fg.endNodes,startNodes,tM,lenObs);
        plhs[0] = getMexArray(maxPath);
        plhs[1] = getMexArray(fg.strainParents);
    }else{
        std::vector<std::vector<int> >maxPath;
        plhs[0] = getMexArray(maxPath);
        plhs[1] = getMexArray(fg.strainParents);
    }

    return;
}
