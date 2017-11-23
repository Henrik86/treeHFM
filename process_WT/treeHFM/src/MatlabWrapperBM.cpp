//
//  MatlabWrapperBM.cpp
//  FGBMA
//
//  Created by Henrik Failmezger on 16.01.15.
//  Copyright (c) 2015 Henrik Failmezger. All rights reserved.
//

#include <stdio.h>
#include "BaumWelch.h"

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
    mxArray * mx = mxCreateDoubleMatrix(v.size(),v.at(0).size(), mxREAL);
    for (int i=0;i<v.size();i++){
        for (int j=0;j<v.at(i).size();j++){
            //std::copy(v.at(i).begin(), v.at(i).end(), arr);
            mxGetPr(mx)[i+(j*v.size())]=v.at(i).at(j);
        }
    }
    return mx;
}
mxArray * getMexArray(std::vector<std::vector<std::vector<double> > >v,int dim1, int dim2, int dim3){
    mwSize dimS [3];
    dimS[0]=dim1;
    dimS[1]=dim2;
    dimS[2]=dim3;
    mxArray * mx = mxCreateNumericArray(3,dimS,mxDOUBLE_CLASS, mxREAL);
    int counter=0;
    for (int k=0;k<v.size();k++){
        for (int i=0;i<v.at(k).size();i++){
            for (int j=0;j<v.at(k).at(i).size();j++){
                //std::copy(v.at(i).begin(), v.at(i).end(), arr);
                mxGetPr(mx)[counter]=v.at(k).at(i).at(j);
                counter++;
            }
        }
    }
    //
    /*
     for(int k=0;k<mxGetDimensions(prhs[11])[2];k++){
     for(int i=0;i<mxGetDimensions(prhs[11])[0];i++){
     for(int j=0;j<mxGetDimensions(prhs[11])[1];j++){
     SigmaInit.at(k).at(i).at(j)=pData11[i+j*mxGetDimensions(prhs[11])[0]+k*mxGetDimensions(prhs[11])[0]* mxGetDimensions(prhs[11])[1]];
     }
     }
     }
     */
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

//void mexFunction( int nlhs, mxArray *plhs[],
//                  int nrhs, const mxArray *prhs[] ){
void mexFunction( int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[] ){
    //DEFINITION////////
    /*
     std::vector<std::vector<double> > startProbN;
     std::vector <double> s1=make_vector<double>() << 0.5 << 0.4 << 0.1;
     startProbN.push_back(s1);
     
     
     std::vector<std::vector <double> > transProbSeq;
     std::vector <double> v1=make_vector<double>() << 0.25<<0.5<<0.25;
     std::vector <double> v2=make_vector<double>() << 0.5<<0.25<<0.25;
     std::vector <double> v3=make_vector<double>() << 0.33<<0.33<<0.33;
     transProbSeq.push_back(v1);
     transProbSeq.push_back(v2);
     transProbSeq.push_back(v3);
     
     
     //double stateIndicesSingle[][2]={{1,1},{1,2},{1,3},{2,1},{2,2},{2,3},{3,1},{3,2},{3,3}};
     std::vector<std::vector <int> > stateIndicesSingle;
     std::vector <int> v4=make_vector<int>() << 1 << 1;
     std::vector <int> v5=make_vector<int>() << 1 << 2;
     std::vector <int> v6=make_vector<int>() << 1 << 3;
     std::vector <int> v7=make_vector<int>() << 2 << 1;
     std::vector <int> v8=make_vector<int>() << 2 << 2;
     std::vector <int> v9=make_vector<int>() << 2 << 3;
     std::vector <int> v10=make_vector<int>() << 3 << 1;
     std::vector <int> v11=make_vector<int>() << 3 << 2;
     std::vector <int> v12=make_vector<int>() << 3 << 3;
     stateIndicesSingle.push_back(v4);
     stateIndicesSingle.push_back(v5);
     stateIndicesSingle.push_back(v6);
     stateIndicesSingle.push_back(v7);
     stateIndicesSingle.push_back(v8);
     stateIndicesSingle.push_back(v9);
     stateIndicesSingle.push_back(v10);
     stateIndicesSingle.push_back(v11);
     stateIndicesSingle.push_back(v12);
     
     
     
     
     //int summationIndices[][3]={{1,2,3},{4,5,6},{7,8,9}};
     std::vector<std::vector <int> > summationIndices;
     std::vector <int> v13=make_vector<int>() << 1 << 2 << 3;
     std::vector <int> v14=make_vector<int>() << 4 << 5 << 6;
     std::vector <int> v15=make_vector<int>() << 7 << 8 << 9;
     summationIndices.push_back(v13);
     summationIndices.push_back(v14);
     summationIndices.push_back(v15);
     summationIndices=subEl(summationIndices,1);
     //
     std::vector<std::vector <double> > transProbDiv;
     std::vector <double> v16=make_vector<double>() << 0.0625<<0.0625<<0.0625<<0.0625<<0.0625<<0.0625<<0.0625<<0.0625<<0.5;
     std::vector <double> v17=make_vector<double>() << 0.111<<0.111<<0.111<<0.111<<0.111<<0.111<<0.111<<0.111<<0.111;
     std::vector <double> v18=make_vector<double>() << 0.111<<0.111<<0.111<<0.111<<0.111<<0.111<<0.111<<0.111<<0.111;
     transProbDiv.push_back(v16);
     transProbDiv.push_back(v17);
     transProbDiv.push_back(v18);
     //
     
     //
     std::vector<std::vector <double> > emProbN;
     std::vector <double> v19=make_vector<double>() << 0.8<<0.1<<0.1;
     std::vector <double> v20=make_vector<double>() << 0.1<<0.8<<0.1;
     std::vector <double> v21=make_vector<double>() << 0.1<<0.1<<0.8;
     emProbN.push_back(v19);
     emProbN.push_back(v20);
     emProbN.push_back(v21);
     //
     
     
     int nHStates=3;
     int nDStates=9;
     
     int nEStates=3;
     
     int stateRep=3;
     
     //1  2  3  4  5  6  7  8  9
     //rr rg rb gr gg gb br bg bb
     //rr gr br rg gg bg rb gb bb
     //int indicesXD[]={1,4,7,2,5,8,3,6,9};
     std::vector <int> indicesXD=make_vector<int>() << 1 << 4 << 7 << 2 << 5 << 8 << 3 << 6 << 9;
     indicesXD= subEl(indicesXD,1); //bring them to c++ indices
     
     std::vector <int> observationSequences=make_vector<int>() <<1<<2<<1<<1<<1<<2;
     std::vector <int> nodeIndices=make_vector<int>()<< 1<< 2<< 3<<4<<5<<6;
     std::vector <int> parents=make_vector<int>() <<0<<1<< 2<<2<< 3<<4;
     
     
     //
     
     
     std::vector<std::vector<double> > transProbDiv;
     std::vector<std::vector<double> > emProbN;
     std::vector<int> observationSequences;
     std::vector<int> parents;
     int vSize;
     std::vector <int> indicesXD;
     
     // std::vector<std::vector<double> > startProbN;
     */
    // startProbN,transProbSeqN,transProbDivUD,probO,observationSequences{os},nodeIndices{os}(:,2),nodeIndices{os}(:,3),summationIndices,indicesXD,'continous'
    //////////////////////// TYPE //////////////////////////////////////////////
    mxChar * pData9 = mxGetChars(prhs[9]);
    char type = pData9[0];
    ////////////////////////////////////////////////////////
    clock_t startLoadData, finishLoadData;
    startLoadData=clock();
    
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
    
    ////////////observation sequences//////////////////////////////////////
    mwSize nDims4 = mxGetNumberOfDimensions(prhs[4]);
    double * pData4 = mxGetPr(prhs[4]);
    mwSize total_num_of_cells;
    total_num_of_cells = mxGetDimensions(prhs[4])[1];
    std::vector <observationSequence> allObservationSequences(total_num_of_cells);
    std::vector <observationSequenceDiscrete> allObservationSequencesDiscrete(total_num_of_cells);
    const mxArray *cell_element_ptr;
    int dimensionData=0;
    if(type=='c'){
    for(int k=0;k<total_num_of_cells;k++){
        cell_element_ptr = mxGetCell(prhs[4], k);
        double * cellData= mxGetPr(cell_element_ptr);
        mwSize cellDataDims = mxGetNumberOfDimensions(cell_element_ptr);
        //
        
        //
        observationSequence obsSequence;
        std::vector<std::vector <double> >  observationMatrix(mxGetDimensions(cell_element_ptr)[0],mxGetDimensions(cell_element_ptr)[1]);
        //
        for(int i=0;i<mxGetDimensions(cell_element_ptr)[0];i++){
            for(int j=0;j<mxGetDimensions(cell_element_ptr)[1];j++){
                observationMatrix.at(i).at(j)=cellData[i+j*mxGetDimensions(cell_element_ptr)[0]];
            }
        }
        dimensionData=observationMatrix.size();
        obsSequence.data=observationMatrix;
        allObservationSequences.at(k)=obsSequence;

    }
    }else{
        for(int k=0;k<total_num_of_cells;k++){
            cell_element_ptr = mxGetCell(prhs[4], k);
            double * cellData= mxGetPr(cell_element_ptr);
            mwSize cellDataDims = mxGetNumberOfDimensions(cell_element_ptr);
        //
        //
            observationSequenceDiscrete obsSequence;
            std::vector <int>  observationMatrix(mxGetDimensions(cell_element_ptr)[1]);
        //
            for(int i=0;i<mxGetDimensions(cell_element_ptr)[1];i++){
                observationMatrix.at(i)=cellData[i];
            }
            dimensionData=observationMatrix.size();
            obsSequence.data=observationMatrix;
           
            allObservationSequencesDiscrete.at(k)=obsSequence;
        }
    }
    ///////////////////////////////////////////////////////////////////////
    
    ////////////node indices//////////////////////////////////////
    mwSize nDims5 = mxGetNumberOfDimensions(prhs[5]);
    double * pData5 = mxGetPr(prhs[5]);
    mwSize total_num_of_cellsN;
    total_num_of_cellsN = mxGetDimensions(prhs[5])[1];
    std::vector <nodeIndices> allNodeIndices(mxGetDimensions(prhs[5])[1]);
    const mxArray *cell_element_ptrN;
    for(int i=0;i<total_num_of_cellsN;i++){
        cell_element_ptrN = mxGetCell(prhs[5], i);
        double * cellData= mxGetPr(cell_element_ptrN);
        mwSize cellDataDims = mxGetNumberOfDimensions(cell_element_ptrN);
        //
        //
        nodeIndices obsNodeIndices;
        std::vector <int>   observationNodeIndices(mxGetDimensions(cell_element_ptrN)[0]);
        //
        for(int j=0;j<mxGetDimensions(cell_element_ptrN)[0];j++){
            observationNodeIndices.at(j)=cellData[j];
        }
        obsNodeIndices.data=observationNodeIndices;
        allNodeIndices.at(i)=obsNodeIndices;
    }
    ////////////parents//////////////////////////////////////
    mwSize nDims6 = mxGetNumberOfDimensions(prhs[6]);
    double * pData6 = mxGetPr(prhs[6]);
    mwSize total_num_of_cellsP;
    total_num_of_cellsP = mxGetDimensions(prhs[6])[1];
    std::vector <parentIndices> allParents(mxGetDimensions(prhs[6])[1]);
    const mxArray *cell_element_ptrP;
    for(int i=0;i<total_num_of_cellsP;i++){
        cell_element_ptrP = mxGetCell(prhs[6], i);
        double * cellData= mxGetPr(cell_element_ptrP);
        mwSize cellDataDims = mxGetNumberOfDimensions(cell_element_ptrP);
        //
        //
        parentIndices obsParent;
        std::vector <int>   observationParent(mxGetDimensions(cell_element_ptrP)[0]);
        //
        for(int j=0;j<mxGetDimensions(cell_element_ptrP)[0];j++){
            observationParent.at(j)=cellData[j];
        }
        obsParent.data=observationParent;
        allParents.at(i)=obsParent;
    }
    ///////////////////////////////////////////////////////////////////////
    
    ///////////////////////////////////////////////////////////////////////
    ////////////summation indices//////////////////////////////////////
    mwSize nDims7 = mxGetNumberOfDimensions(prhs[7]);
    double * pData7 = mxGetPr(prhs[7]);
    std::vector <std::vector <int> > summationIndices(mxGetDimensions(prhs[7])[0], std::vector<int>(mxGetDimensions(prhs[7])[1]));
    for(int i=0;i<mxGetDimensions(prhs[7])[0];i++){
        for(int j=0;j<mxGetDimensions(prhs[7])[1];j++){
            summationIndices.at(i).at(j)=pData7[i+j*mxGetDimensions(prhs[7])[0]];
        }
    }
    summationIndices=subEl(summationIndices,1);
    ///////////////////////////////////////////////////////////////////////
    ////////////indices XD//////////////////////////////////////
    mwSize nDims8 = mxGetNumberOfDimensions(prhs[8]);
    double * pData8 = mxGetPr(prhs[8]);
    std::vector <int> indicesXD(mxGetDimensions(prhs[8])[1]);
    for(int i=0;i<mxGetDimensions(prhs[8])[1];i++){
        indicesXD.at(i)=pData8[i];
    }
    indicesXD= subEl(indicesXD,1); //bring them to c++ indices
    //////////// OBS PROB FOR CONTINOUS ///////////////////////////////////
    
    //////////////////////// StateIndicesSingle //////////////////////////////////////////////
    mwSize nDims10 = mxGetNumberOfDimensions(prhs[10]);
    double * pData10 = mxGetPr(prhs[10]);
    std::vector<std::vector <int> > stateIndicesSingle(mxGetDimensions(prhs[10])[0], std::vector<int>(mxGetDimensions(prhs[10])[1]));
    for(int i=0;i<mxGetDimensions(prhs[10])[0];i++){
        for(int j=0;j<mxGetDimensions(prhs[10])[1];j++){
            stateIndicesSingle.at(i).at(j)=pData10[i+j*mxGetDimensions(prhs[10])[0]];
        }
    }
    finishLoadData=clock();
    //////////////////////// SigmaInit //////////////////////////////////////////////
    mwSize nDims11 = mxGetNumberOfDimensions(prhs[11]);
    double * pData11 = mxGetPr(prhs[11]);
    int sdim1,sdim2,sdim3;
    if(type=='c'){
        sdim1=mxGetDimensions(prhs[11])[2];
        sdim2=mxGetDimensions(prhs[11])[1];
        sdim3=mxGetDimensions(prhs[11])[0];
    }else{
        sdim1=1;
        sdim2=1;
        sdim3=1;
    }
    std::vector<std::vector<std::vector<double> > > SigmaInit(sdim1, std::vector< std::vector<double> >(sdim2, std::vector<double> (sdim3)));
    //
    if(type=='c'){
    //
         for(int k=0;k<mxGetDimensions(prhs[11])[2];k++){
             for(int i=0;i<mxGetDimensions(prhs[11])[0];i++){
                 for(int j=0;j<mxGetDimensions(prhs[11])[1];j++){
                     SigmaInit.at(k).at(i).at(j)=pData11[i+j*mxGetDimensions(prhs[11])[0]+k*mxGetDimensions(prhs[11])[0]* mxGetDimensions(prhs[11])[1]];
                 }
             }
         }
     }
    ///////////////////////////////////////////////////////////////////////
    //////////////////////// MuInit //////////////////////////////////////////////
    mwSize nDims12 = mxGetNumberOfDimensions(prhs[12]);
    double * pData12 = mxGetPr(prhs[12]);
    int mdim1,mdim2;
    if(type=='c'){
            mdim1=mxGetDimensions(prhs[11])[0];
            mdim2=mxGetDimensions(prhs[11])[1];
    }else{
            mdim1=1;
            mdim2=1;
    }
    std::vector<std::vector<double> > muN(mdim1, std::vector<double>(mdim2));
    //
    //
     if(type=='c'){
         for(int i=0;i<mxGetDimensions(prhs[12])[0];i++){
             for(int j=0;j<mxGetDimensions(prhs[12])[1];j++){
                 muN.at(i).at(j)=pData12[i+j*mxGetDimensions(prhs[12])[0]];
             }
         }
     }
    ///////////////////////////////////////////////////////////////////////
    if(type=='c'){
        outputC bmOutput=baumWelch(allObservationSequences,  startProbN, transProbSeqN,  transProbDiv,  emProbN, muN, SigmaInit,  allNodeIndices, summationIndices, stateIndicesSingle, indicesXD,  allParents,  type);
       plhs[0] = getMexArray(bmOutput.transProbSeq);
       plhs[1] = getMexArray(bmOutput.transProbDiv);
       plhs[2] = getMexArray(bmOutput.prior);
       plhs[3] = getMexArray(bmOutput.allLL);
       plhs[4] = getMexArray(bmOutput.allLLSingle);
       plhs[5] = getMexArray(bmOutput.mu);
       plhs[6] = getMexArray(bmOutput.sigma,dimensionData,dimensionData,bmOutput.transProbSeq.size());
    }else{
        outputD bmOutputDiscrete=baumWelchDiscrete(allObservationSequencesDiscrete,  startProbN, transProbSeqN,  transProbDiv,  emProbN, allNodeIndices, summationIndices, stateIndicesSingle, indicesXD,  allParents,  type);
        
        plhs[0] = getMexArray(bmOutputDiscrete.transProbSeq);
        plhs[1] = getMexArray(bmOutputDiscrete.transProbDiv);
        plhs[2] = getMexArray(bmOutputDiscrete.prior);
        plhs[3] = getMexArray(bmOutputDiscrete.emProb);
        plhs[4] = getMexArray(bmOutputDiscrete.allLL);
        plhs[5] = getMexArray(bmOutputDiscrete.allLLSingle);

    }
    ///////////////////////////////////////////////////////////////////////
    
    
    return;
}
