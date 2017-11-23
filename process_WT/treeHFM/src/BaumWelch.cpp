//============================================================================
// Name        : FactorGraph.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================



#include "BaumWelch.h"

using namespace std;


//FUNCTIONS:
//print array
//print vector
//print2D array
std::vector<std::vector<std::vector<double> > > calculateEijSeq(factorGraph &fgS,std::vector<std::vector<double> >  &transProbSeq,std::vector<std::vector<double> > &observationProbabilities,int numberNodes ){
    
    std::vector<std::vector<std::vector<double> > > EijSeq(numberNodes+1);
    int nHStates=transProbSeq.at(0).size();
    //EijDiv=initialise(transProbDivN.at(0).size());
    for (int t=1; t<=numberNodes; t++) {
            std::vector<int> f_neighbors=fgS.variableNodes.at(t).f_neighbors;
            //
            Max  max_f_neighbors=max(f_neighbors);
            int f_neighbor_max=max_f_neighbors.value;
            std::vector<double> msg_down;
            std::vector<double> msg_up;

            if(!fgS.factorNodes.at(f_neighbor_max).division && intersect(fgS.endNodes,t)==0){
                std::vector<std::vector<double> > EijSeqNode(nHStates,std::vector<double> (nHStates));
                std::vector<int> f_indexUp=fgS.factorNodes.at(f_neighbor_max).f_children;
                Max  max_v_neighbors=max(fgS.variableNodes.at(t).v_neighbors);
                int v_neighbors=max_v_neighbors.value;
                //
                if(intersect(fgS.endNodes,v_neighbors)==1){
                    msg_down=fgS.factorNodes.at(f_neighbors.at(0)).msg_out_down.at(findMessage(fgS.factorNodes.at(f_neighbors.at(0)).msg_out_down,t)).msg;
                    msg_up=ones(msg_down.size());
                }else{
                    msg_up=fgS.factorNodes.at(f_indexUp.at(0)).msg_out_up.at(findMessage(fgS.factorNodes.at(f_indexUp.at(0)).msg_out_up,v_neighbors)).msg;
                    msg_down=fgS.factorNodes.at(f_neighbors.at(0)).msg_out_down.at(findMessage(fgS.factorNodes.at(f_neighbors.at(0)).msg_out_down,t)).msg;
                }

                double factor=0.0;
                for(int n1=0;n1<nHStates;n1++){
                     for(int d1=0;d1<nHStates;d1++){
                        factor=factor + msg_down.at(n1) * msg_up.at(d1) * transProbSeq.at(n1).at(d1) * observationProbabilities.at(v_neighbors-1).at(d1);

                    }
                 }
                 for(int n1=0;n1<nHStates;n1++){
                     for(int d1=0;d1<nHStates;d1++){
                        EijSeqNode.at(n1).at(d1)=((msg_down.at(n1)*msg_up.at(d1)*transProbSeq.at(n1).at(d1)*observationProbabilities.at(v_neighbors-1).at(d1))/factor);

                      }
                }

                EijSeq.at(t)=EijSeqNode;
            }else{
                std::vector<std::vector<double> > EijSeqT;
                EijSeq.at(t)=EijSeqT;
            }
    }
    return EijSeq;
}

std::vector<std::vector<std::vector<double> > >  calculateEijDiv(factorGraph &fgS,std::vector<std::vector<double> >  &transProbDivN,std::vector<std::vector<double> > &observationProbabilities,std::vector<std::vector<int> > stateIndicesSingle, int numberNodes, int nHStates, int nDStates){
    std::vector<std::vector<std::vector<double> > > EijDiv;
    
    
    //EijDiv=initialise(transProbDivN.at(0).size());
    for (int t=1; t<=numberNodes; t++) {
            std::vector<int> f_neighbors=fgS.variableNodes.at(t).f_neighbors;
            //
            Max max_f_neighbors=max(f_neighbors);
            int f_neighbors_max=max_f_neighbors.value;
            //
            if(fgS.factorNodes.at(f_neighbors_max).division){
                if(!intersect(fgS.endNodes,t)){
                    //
                    std::vector<std::vector<double> > EijDivNode(nHStates,std::vector<double> (nDStates));
                    //
                    int nHStates=transProbDivN.size();
                    int nDStates=transProbDivN.at(0).size();
                    //
                    Min min_f_neighbors=min(f_neighbors);
                    int f_neighbors_min=min_f_neighbors.value;
                    //
                    std::vector<int>  f_indexUp=fgS.factorNodes.at(f_neighbors_max).f_children;
                    //
                    Min  min_f_indexUp=min(f_indexUp);
                    int f_indexUpPlus=min_f_indexUp.value;
                    //
                    Max  max_f_neighbors=max(f_indexUp);
                    int f_indexUpMinus=max_f_neighbors.value;
                    //
                    std::vector<int> v_neighbors=fgS.factorNodes.at(f_neighbors_max).v_neighbors;
                    //
                    Min min_v_neighbors=min(fgS.variableNodes.at(t).v_neighbors);
                    int v_neighbors_min=min_v_neighbors.value;
                    int v_neighbors_minIndex=min_v_neighbors.index;
                    //
                    std::vector<int> v_neighborsDown=v_neighbors;
                    v_neighbors.erase(v_neighbors.begin()+v_neighbors_minIndex); //delete min v neighbor (neighbor downwards)
                    //
                    Min min_v_neighbors2=min(v_neighbors);
                    int v_neighborsUpPlus=min_v_neighbors2.value;//variable node +
                    //
                    Max  max_v_neighbors=max(v_neighbors);
                    int v_neighborsUpMinus=max_v_neighbors.value;//variable node -
                    //
                    //
                    std::vector<double> msg_upPlus;
                    //
                    std::vector<double> msg_down=fgS.factorNodes.at(f_neighbors_min).msg_out_down.at(findMessage(fgS.factorNodes.at(f_neighbors_min).msg_out_down,t)).msg;
                    if(intersect(fgS.endNodes,v_neighborsUpPlus)){
                        msg_upPlus=ones(msg_down.size());
                    }else{
                        msg_upPlus=fgS.factorNodes.at(f_indexUpPlus).msg_out_up.at(findMessage(fgS.factorNodes.at(f_indexUpPlus).msg_out_up,v_neighborsUpPlus)).msg;
                    }
                    //
                    std::vector<double> msg_upMinus;
                    //
                    if(intersect(fgS.endNodes,v_neighborsUpMinus)){
                        msg_upMinus=ones(msg_down.size());
                    }else{
                        msg_upMinus=fgS.factorNodes.at(f_indexUpMinus).msg_out_up.at(findMessage(fgS.factorNodes.at(f_indexUpMinus).msg_out_up,v_neighborsUpMinus)).msg;
                    }
                    
                    double factor=0.0;
                    for(int n1=0;n1<nHStates;n1++){
                        for(int d1=0;d1<nDStates;d1++){
                            std::vector <int> si=stateIndicesSingle.at(d1);
                            factor=factor+msg_down.at(n1)*msg_upPlus.at(si.at(1-1)-1)*msg_upMinus.at(si.at(2-1)-1)*transProbDivN.at(n1).at(d1)*observationProbabilities.at(v_neighborsUpPlus-1).at(si.at(1-1)-1)*observationProbabilities.at(v_neighborsUpMinus-1).at(si.at(2-1)-1);
                        }
                    }
                    for(int n1=0;n1<nHStates;n1++){
                        for(int d1=0;d1<nDStates;d1++){
                            std::vector <int>si=stateIndicesSingle.at(d1);
                            EijDivNode.at(n1).at(d1)=((msg_down.at(n1)*msg_upPlus.at(si.at(1-1)-1)*msg_upMinus.at(si.at(2-1)-1)*transProbDivN.at(n1).at(d1)*observationProbabilities.at(v_neighborsUpPlus-1).at(si.at(1-1)-1)*observationProbabilities.at(v_neighborsUpMinus-1).at(si.at(2-1)-1))/factor);
                        }
                    }
                    EijDiv.push_back(EijDivNode);
            }
        }
    
    }
    return EijDiv;
}

std::vector<std::vector<double> > calculateYpsilon(factorGraph &fgS,std::vector<std::vector<double> > &observationProbabilities,int numberNodes){
    std::vector<std::vector<double> > ypsilonAllProb(numberNodes+1,std::vector<double>(fgS.factorNodes.at(0).msg_out_down.size()));
   
    for (int t=1; t<=numberNodes; t++) {
            std::vector<double> Ypsilon;
            double factor=0.0;
            std::vector<int> f_neighbors=fgS.variableNodes.at(t).f_neighbors;
            std::vector<double> msg_down;
            std::vector<double> msg_up;
            if(intersect(fgS.endNodes,t)){
             msg_down=fgS.factorNodes.at(f_neighbors.at(0)).msg_out_down.at(findMessage(fgS.factorNodes.at(f_neighbors.at(0)).msg_out_down,t)).msg;
             msg_up=ones(msg_down.size());
            }else{
             msg_up=fgS.factorNodes.at(f_neighbors.at(1)).msg_out_up.at(findMessage(fgS.factorNodes.at(f_neighbors.at(1)).msg_out_up,t)).msg;
             msg_down=fgS.factorNodes.at(f_neighbors.at(0)).msg_out_down.at(findMessage(fgS.factorNodes.at(f_neighbors.at(0)).msg_out_down,t)).msg;
            }
            double sumMSG=sum(multEl(msg_down,msg_up));
            for(int n1=0;n1<msg_down.size();n1++){
                Ypsilon.push_back(msg_down.at(n1)*msg_up.at(n1)/sumMSG);
            }
        
            ypsilonAllProb.at(t-1)=Ypsilon;
    }
    return ypsilonAllProb;
}

std::vector<std::vector<std::vector<double> > >  calculateSigma(std::vector<observationSequence> allObservationSequences, int nHStates, std::vector<ypsilon> allYpsilons,std::vector <int> skipSequences,std::vector<std::vector<double> > muN,double covarianceFaktor, int dimData,std::vector<std::vector<std::vector<double> > >  covOld ){
    //
    std::vector<std::vector<std::vector<double> > > SigmaN(nHStates) ;
    int nObs =allObservationSequences.size();
    //nObs=1;
    //
   
    for(int n1=0;n1<nHStates;n1++){
        double w=0; //%sum(ypsilon)
        std::vector<std::vector<double> >  m;
        std::vector<double> muState=getColumn(muN,n1);
        bool mDef=0;
        for(int os=0;os<nObs;os++){
            
            if(allObservationSequences.at(os).data.size()>0 && skipSequences.at(os)==0){
                 std::vector<std::vector<double> >ypsilon=allYpsilons.at(os).data;
                 std::vector<std::vector<double> > dataS=allObservationSequences.at(os).data;

                 std::vector<std::vector<double> > dataComponent=dataS;
                 std::vector<double> ypsilonState=getColumn(ypsilon,n1);
                if(n1==0){
                 //printVector(ypsilonState);
                }
                 std::vector<std::vector<double> > ypsilonComponent=repmat(ypsilonState,dimData,1);
                 //
                 w=w+sum(ypsilonState); //%FORMEL NACHPRUEFEN !!!
                 if(mDef==0){//first observation
                     //printVector((multEl(ypsilonComponent,dataComponent)));
                     //printVector(matrixMultiplication(multEl(ypsilonComponent,dataComponent),transpose(dataComponent)));
                     m=matrixMultiplication(multEl(ypsilonComponent,dataComponent),transpose(dataComponent));
                     mDef=1;
                 }else{
                     m=addEl(m,matrixMultiplication(multEl(ypsilonComponent,dataComponent),transpose(dataComponent)));
                 }
               
                
                //printVector(multEl(multEl(ypsilonComponent,dataComponent),dataComponent));
            }
        }
        if(w==0){//%%%Is this right?????
            w=1;
        }
        std::vector<std::vector<double> > matrixTemp;
        matrixTemp.push_back(muState);

        SigmaN.at(n1)=subEl(divEl(m,w),matrixMultiplication(transpose(muState),matrixTemp));
                      //
    }
    covarianceFaktor=0.08;
    SigmaN=addEl(SigmaN,covarianceFaktor);
    SigmaN=addEl(multEl(SigmaN,0.8),multEl(covOld,0.2));
    return SigmaN;
}
std::vector<double> calculatePI(int nHStates,std::vector<ypsilon> allYpsilon, std::vector <int> skipSequences){
    std::vector<double> piK(nHStates);
    int nObs=skipSequences.size();
    //nObs=1;
    for (int n=0; n<nHStates;n++){
        piK.at(n)=0;
        for(int os=0;os<nObs;os++){
            if(skipSequences.at(os)==0){
                      piK.at(n)=piK.at(n)+allYpsilon.at(os).data.at(0).at(n); //%./sum(allYpsilon{os}(1,:))
            }
        }
        piK.at(n)=piK.at(n);//%/nHStates;
    }
    piK=divEl(piK,sum(piK));
    return piK;
}
std::vector< std::vector<double> > calculateBJ(std::vector<observationSequenceDiscrete> allObservationSequences,int nHStates,std::vector<ypsilon> allYpsilon, std::vector <int> skipSequences, std::vector<int> oStates){
    
    //
    int nOStates= oStates.size();
    std::vector< std::vector<double> > bJ(nHStates,std::vector<double>(nOStates));
    int nObs=allObservationSequences.size();
    //nObs=1;
    std::vector<double> factor(nHStates);
    for (int j=0; j<nHStates;j++){
        for(int os=0;os<nObs;os++){
            for(int t=0; t<allObservationSequences.at(os).data.size(); t++){
                factor.at(j)=factor.at(j)+allYpsilon.at(os).data.at(t).at(j);
            }
        }
    }
    std::vector< std::vector<double> > numerator(nOStates,std::vector<double>(nHStates));
    for (int j=0; j<nHStates;j++){
        for (int k=0; k<nOStates;k++){
            for(int os=0;os<nObs;os++){
                for(int t=0; t<allObservationSequences.at(os).data.size(); t++){
                    if(allObservationSequences.at(os).data.at(t) == oStates.at(k)){
                        numerator.at(k).at(j)=numerator.at(k).at(j)+allYpsilon.at(os).data.at(t).at(j);
                    }
                }
            }
        }
    }
    for (int j=0; j<nHStates;j++){
        for (int k=0; k<nOStates;k++){
            bJ.at(j).at(k)=numerator.at(k).at(j)/factor.at(j);
        }
    }
    return bJ;
}
std::vector<std::vector<double> > calculateTransProbSeq(int nHStates, std::vector <int> skipSequences,std::vector <eijSeq>  allEijSeq ){
    //%%%% estimate the new transition matrix for sequences %%%%%
    std::vector<std::vector<double> > transProbSeqN_tmp=zeros(nHStates,nHStates);
    //%
    int nObs=skipSequences.size();
    //nObs=1;
    for(int os=0; os<nObs;os++){
        if( skipSequences.at(os)==0){
            std::vector< std::vector<std::vector<double> > > Eij=allEijSeq.at(os).data;
            for(int n1=0;n1<nHStates;n1++){
                for (int d1=0;d1<nHStates;d1++){
                    double numer=0;
                    for(int t=0; t<Eij.size(); t++){
                        if(Eij.at(t).size()>0){
                            numer=numer+Eij.at(t).at(n1).at(d1);
                            if(n1==2 && d1==0){
                            }
                        }
                    }
                    
                    transProbSeqN_tmp.at(n1).at(d1)=transProbSeqN_tmp.at(n1).at(d1)+numer;
                }
            }
        }
    }
   std::vector<double> ZS = sum(transProbSeqN_tmp,1);
   std::vector<double>  SS = addEl(ZS,logicalMatrix(ZS,0));//%%replace 0 by 1
   std::vector<std::vector<double> >  normS = repmat(SS,transProbSeqN_tmp.at(0).size(),2);
   std::vector<std::vector<double> > transProbSeqN=divEl(transProbSeqN_tmp,normS);
   return transProbSeqN;
}
std::vector<std::vector<double> > calculateTransProbDiv(int nHStates,int nDStates, std::vector <int> skipSequences,std::vector <eijDiv>  allEijDiv){
    std::vector<std::vector<double> > transProbDivN_tmp=zeros(nHStates,nDStates);
    //%
    int nObs=skipSequences.size();
    //nObs=1;
    //
    if(allEijDiv.size()>0){
        for(int os=0;os < nObs;os++)
            if( skipSequences.at(os)==0){
                std::vector< std::vector<std::vector<double> > >  Eij=allEijDiv.at(os).data;
                for(int n1=0;n1< nHStates;n1++){
                    for(int d1=0;d1<nDStates;d1++){
                            double numer=0;
                            double denom=0;
                            for(int t=0;t< Eij.size();t++){
                              numer=numer+Eij.at(t).at(n1).at(d1);
                            }
                            transProbDivN_tmp.at(n1).at(d1)=transProbDivN_tmp.at(n1).at(d1)+numer;
                    }
                }
            }
    }
    std::vector<double> ZD = sum(transProbDivN_tmp,1);
    std::vector<double> SD = addEl(ZD ,logicalMatrix(ZD,0));//%%replace 0 by 1
    std::vector<std::vector<double> > normD = repmat(SD,transProbDivN_tmp.at(0).size(),2);
    std::vector<std::vector<double> >  transProbDivN=divEl(transProbDivN_tmp,normD);
    return transProbDivN;
}
std::vector<std::vector<double> > calculateMu(std::vector <observationSequence> allObservationSequences,int nHStates, std::vector <int> skipSequences,std::vector <ypsilon> allYpsilon, int dimData){
    std::vector<std::vector<double> > muN(dimData,std::vector<double>(nHStates));
    
    int nObs=allObservationSequences.size();
    //int nObs=1;
    for (int s=0;s<nHStates; s++){
        std::vector<double> numerator =zeros(dimData);
        double factor=0;
        for (int os=0;os<nObs; os++){
             if(allObservationSequences.at(os).data.size()>0 && skipSequences.at(os)==0){
                std::vector<std::vector<double> >  ypsilon=allYpsilon.at(os).data;
                std::vector<double>  ypsilonState=getColumn(ypsilon,s);
                std::vector<std::vector<double> >  ypsilonComponent=repmat(ypsilonState,dimData,2);
                 //
                std::vector<std::vector<double> > data=allObservationSequences.at(os).data;
                data=transpose(data);
                multEl(ypsilonComponent,data);
                std::vector<double> ss=sum(multEl(ypsilonComponent,data),2);
                addEl(numerator,sum(multEl(ypsilonComponent,data),1));
                numerator=addEl(numerator,sum(multEl(ypsilonComponent,data),2));
                factor=factor+sum(ypsilonState);
             }
         }
                             ///%muN(:,s)=(numerator/factor)';
                             //
        if (factor==0){ //%%%Is this right?????
            factor=1;
        }
        muN=replaceColumn(muN,divEl(numerator,factor),s);
                             
    }
    return muN;
}
outputC baumWelch( std::vector <observationSequence> allObservationSequences, std::vector<std::vector<double> > startProbN, std::vector<std::vector<double> > transProbSeqN, std::vector<std::vector<double> > transProbDivN, std::vector<std::vector<double> > emProbN, std::vector<std::vector<double> > muN,  std::vector<std::vector<std::vector<double> > > SigmaN, std::vector<nodeIndices> allNodeIndices, std::vector<std::vector<int> > summationIndices, std::vector<std::vector<int> > stateIndicesSingle, std::vector<int> indicesXD, std::vector<parentIndices > allParentIndices, int type){
    

//'BAUM WELCH CMO';
std::vector<double> allLL;
std::vector<std::vector<double> > allLLSingle;
int converged=0;
double ll;
double ll_old;
int iterations=0;
int nHStates=transProbSeqN.size();
int nDStates=transProbDivN.at(0).size();
int vSize=transProbSeqN.size();
int counterIt=0;
int counter=0;
double eps=2.2204e-16;
std::vector<double> probAllRuns;
//%
int nObs=allObservationSequences.size();
//%
std::vector<ypsilon> allYpsilon(nObs);
std::vector<eijSeq> allEijSeq(nObs);
std::vector<eijDiv> allEijDiv(nObs);
//
std::vector<int> skipSequences=zerosINT(allObservationSequences.size());//%if 1 skip this sequence (division in 3 children)
double covarianceFaktor=0.1;
while(!converged){
    
    int nStart=1;
    //%%%%SCALING NOCHMAL UEBERDENKEN%%%%%%%
    
    int osCounter=0;
    int dimensionData;
    //nObs=1;
    std::vector<double> allProbObs(nObs,NAN);
    for(int osCounter=0;osCounter<nObs;osCounter++){
        std::vector<double>  tPP=tabulate(allParentIndices.at(osCounter).data);
        Max maxParentsInd;
        maxParentsInd = max(tPP);
        //
        int earlyDivision=0;
        if(maxParentsInd.value >2 ||  countIndex(allParentIndices.at(osCounter).data,1)>1  ){
            earlyDivision=1;
            skipSequences.at(osCounter)=1;
        }
        if(allObservationSequences.at(osCounter).data.size()==0){
            skipSequences.at(osCounter)=1;
        }
        if(allObservationSequences.at(osCounter).data.size()>0 && earlyDivision==0){

                 dimensionData=SigmaN.at(0).size();
                 int* start_d = (int*)malloc(sizeof(int)*dimensionData);
                 int o;
                 for (o = 0; o < dimensionData; o++) {
                     start_d[o] = o;
                 }
                  
                  std::vector<std::vector<double> > probO(allObservationSequences.at(osCounter).data.at(0).size(),std::vector<double>(nHStates));
                    //int s=0;
                    //ParamContainerEmissions* multGParams= new ParamContainerEmissions(vectorToArray(getColumn(muN,s)), vectorToArray(SigmaN.at(s)), 0, dimensionData, start_d, 1, 1);
                  //MultivariateGaussian mG=MultivariateGaussian(multGParams);
            
                    for(int d=0;d< allObservationSequences.at(osCounter).data.at(0).size(); d++){
                      for(int s=0; s < nHStates; s++){
                          ParamContainerEmissions* multGParams=new ParamContainerEmissions(vectorToArray(getColumn(muN,s)), vectorToArray(SigmaN.at(s)), 0, dimensionData, start_d, 1, 1);
                          MultivariateGaussian* mG= new MultivariateGaussian(multGParams);
                          double prob = mG->calcEmissionProbability(vectorToArray1D(getColumn(allObservationSequences.at(osCounter).data,d)), 0, 0);
                          probO.at(d).at(s)=prob;
                          delete mG;
                        }
                    }
                  free (start_d);
                  int lenObs=allObservationSequences.at(osCounter).data.at(0).size();
                  int obsProbSize=probO.size();
                   /*
                  dimData=nrows(observationSequences{os});
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  %mex CPP/FactorGraphOpt.cpp CPP/helpFunctions.cpp
                  */
            
                factorGraph fg;
                createHMT(fg,allNodeIndices.at(osCounter).data, allParentIndices.at(osCounter).data,summationIndices,0,vSize,indicesXD,type,probO,transProbSeqN.size(),transProbDivN.at(0).size());
                std::vector<int>endNodes=fg.endNodes;
                std::vector<int> startNodes;
                startNodes.push_back(1);
        
                //
                std::vector< std::vector<std::vector <double> > > transProbDivUD;
                transProbDivUD.push_back(transProbDivN);
                transProbDivUD.push_back(arrangeItems(transProbDivN,indicesXD));
                transitionMatrices tM;
                tM.probDiv=transProbDivUD;
                std::vector< std::vector<std::vector<double> > > startProbT;
                startProbT.push_back(startProbN);
                tM.probStart=startProbT;
                std::vector< std::vector<std::vector<double> > > transProbSeqT;
                transProbSeqT.push_back(transProbSeqN);
                tM.probSeq=transProbSeqT;
                //
                factorGraph fgN=sumProductAlgorithm(fg,fg.endNodes,startNodes,tM);
                 /*
                  if(~isempty(strainParents))
                  tP=tabulate(strainParents(:,1));
                  tP=tP(:,2);
                  else
                  tP=1;
                  end
        */
                  int numberNodes=fgN.factorNodes.size();
                  std::vector<double> tP;
                  if(fgN.strainParents.size()>0){
                     tP=tabulate(getColumn(fgN.strainParents,0));
                  }else{
                      tP.push_back(0);
                  }

                  int endNodeCheck=0;
                  for(int e=0;e<endNodes.size();e++){
                      int eN=endNodes.at(e);
                  }
            
                  //emptyCells = cellfun('isempty', EijSeq); //%delete the empty cells (multiple observations)
                  //EijSeqT(emptyCells) = [];
                  //%%%%%%%%%%%%%%%%%%%%%%%%%
                  //~isempty(find(isnan(fnN(2).msg_out_up{1}.msg))) ||
                  Max maxParents;
                  maxParents = max(tP);
                  if(maxParents.value >=3 ||  endNodeCheck ){
                      skipSequences.at(osCounter)=1;
                      continue;
                  }else{
                      
                      int N_v=allObservationSequences.at(osCounter).data.at(0).size();
                      //emptyCells = cellfun('isempty', EijSeq); //delete the empty cells (multiple observations)
                      //EijSeq(emptyCells) = [];
                      //%%%%%
                     // emptyCells = cellfun('isempty', EijDiv); //delete the empty cells (multiple observations)
                      //
                      eijSeq eijSeqObs;
                      eijSeqObs.data=calculateEijSeq(fgN,transProbSeqN,probO,obsProbSize);
                      allEijSeq.at(osCounter)=eijSeqObs;
                      //
                      eijDiv eijDivObs;
                      eijDivObs.data = calculateEijDiv(fgN,transProbDivN,probO,stateIndicesSingle,obsProbSize,nHStates,nDStates);
                      allEijDiv.at(osCounter)=eijDivObs;
                      //
                      ypsilon ypsilonObs;
                      ypsilonObs.data=calculateYpsilon(fgN,probO,obsProbSize);
                      allYpsilon.at(osCounter)=ypsilonObs;
                 }
         
            /*
           
             */
                allProbObs.at(osCounter)=(log(sum(multEl(fgN.factorNodes.at(1).msg_out_down.at(0).msg,fgN.factorNodes.at(2).msg_out_up.at(0).msg)))+fgN.factorNodes.at(1).acc_logS_msg_out_down.at(0).msg+fgN.factorNodes.at(2).acc_logS_msg_out_up.at(0).msg);
                 // %allProbObs{os}=log(sum(fnO(1).msg_out_down{1}.*fnO(2).msg_out_up{1}))+fnO(1).acc_logS_msg_out_down(1)+fnO(2).acc_logS_msg_out_up(1);
            }
        }
    
        if(counterIt>0){
            ll_old=ll;
        }
        ll=0;
        for(int osCounter=0;osCounter<nObs;osCounter++){
            if(allObservationSequences.at(osCounter).data.size()>0 && skipSequences.at(osCounter)==0){
                ll=ll+allProbObs.at(osCounter);
            }
        }
        allLL.push_back(ll);
        allLLSingle.push_back(allProbObs);
        if(counterIt>0){
            if((abs(ll_old-ll)/(abs(ll_old) + abs(ll) + eps))<1e-04 || counterIt>30){
                converged=1;
                break;
            }
        }
        std::vector <double> piN=calculatePI(nHStates,  allYpsilon, skipSequences);
        startProbN.at(0)=piN;
 
        transProbSeqN=calculateTransProbSeq(nHStates,  skipSequences, allEijSeq );
   
        //

        transProbDivN=calculateTransProbDiv(nHStates,nDStates, skipSequences,allEijDiv);

        muN= calculateMu(allObservationSequences, nHStates, skipSequences,allYpsilon, dimensionData);
        std::vector<std::vector<std::vector<double> > >  SigmaOld= SigmaN;
        SigmaN= calculateSigma(allObservationSequences, nHStates, allYpsilon,skipSequences, muN, covarianceFaktor,dimensionData,  SigmaOld);
        
       
    cout<<"Log-Likelihood: "<<ll<<endl;
    counterIt=counterIt+1;
    iterations=iterations+1;
    }
    outputC bmOutput;
    bmOutput.transProbSeq=transProbSeqN;
    bmOutput.transProbDiv=transProbDivN;
    bmOutput.mu=muN;
    bmOutput.sigma=SigmaN;
    bmOutput.prior=startProbN;
    bmOutput.allLL=allLL;
    bmOutput.allLLSingle=allLLSingle;
    return bmOutput;
}
std::vector<int> findNumberOfDiscreteStates(std::vector <observationSequenceDiscrete> allObservationSequences){
    int maxState=allObservationSequences.at(0).data.at(0);
    for(int i=0; i<allObservationSequences.size(); i++){
            for(int j=0; j<allObservationSequences.at(i).data.size(); j++){
                if(allObservationSequences.at(i).data.at(j)>maxState){
                    maxState=allObservationSequences.at(i).data.at(j);
                }
            }
    }
    std::vector<int> discreteStates=createSequence(maxState+1,1, maxState+1);
    return discreteStates;
}
outputD baumWelchDiscrete( std::vector <observationSequenceDiscrete> allObservationSequences, std::vector<std::vector<double> > startProbN, std::vector<std::vector<double> > transProbSeqN, std::vector<std::vector<double> > transProbDivN, std::vector<std::vector<double> > emProbN,  std::vector<nodeIndices> allNodeIndices, std::vector<std::vector<int> > summationIndices, std::vector<std::vector<int> > stateIndicesSingle, std::vector<int> indicesXD, std::vector<parentIndices > allParentIndices, int type){
   
    std::vector<double> allLL;
    std::vector<std::vector<double> > allLLSingle;
    int converged=0;
    double ll;
    double ll_old;
    int iterations=0;
    int nEStates=emProbN.at(0).size();
    int nHStates=transProbSeqN.size();
    int nDStates=transProbDivN.at(0).size();
    int vSize=transProbSeqN.size();
    int counterIt=0;
    int counter=0;
    double eps=2.2204e-16;
    std::vector<double> probAllRuns;
    //%
    int nObs=allObservationSequences.size();
    std::vector<int>  oStates= findNumberOfDiscreteStates(allObservationSequences);
   
    
    //%
    std::vector<ypsilon> allYpsilon(nObs);
    std::vector<eijSeq> allEijSeq(nObs);
    std::vector<eijDiv> allEijDiv(nObs);
    //
    std::vector<int> skipSequences=zerosINT(allObservationSequences.size());//%if 1 skip this sequence (division in 3 children)
    double covarianceFaktor=0.1;
    while(!converged){
        
        int nStart=1;
        //%%%%SCALING NOCHMAL UEBERDENKEN%%%%%%%
        
        int osCounter=0;
        int dimensionData;
        //nObs=1;
        std::vector<double> allProbObs(nObs,0);
        for(int osCounter=0;osCounter<nObs;osCounter++){
            std::vector<double>  tPP=tabulate(allParentIndices.at(osCounter).data);
            Max maxParentsInd;
            maxParentsInd = max(tPP);
            //
            int earlyDivision=0;
            if(maxParentsInd.value >2 ||  countIndex(allParentIndices.at(osCounter).data,1)>1  ){
                earlyDivision=1;
                skipSequences.at(osCounter)=1;
            }
            if(allObservationSequences.at(osCounter).data.size()>0 && earlyDivision==0){
                
   
                std::vector<std::vector<double> > probO(allObservationSequences.at(osCounter).data.size(),std::vector<double>(nHStates));
                //int s=0;
                //ParamContainerEmissions* multGParams= new ParamContainerEmissions(vectorToArray(getColumn(muN,s)), vectorToArray(SigmaN.at(s)), 0, dimensionData, start_d, 1, 1);
                //MultivariateGaussian mG=MultivariateGaussian(multGParams);
                
                for(int d=0;d< allObservationSequences.at(osCounter).data.size(); d++){
                    for(int s=0; s < nHStates; s++){
                        probO.at(d).at(s)=emProbN.at(s).at(allObservationSequences.at(osCounter).data.at(d)-1);
                    }
                }
                int lenObs=allObservationSequences.at(osCounter).data.size();
                int obsProbSize=probO.size();
                /*
                 dimData=nrows(observationSequences{os});
                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 %mex CPP/FactorGraphOpt.cpp CPP/helpFunctions.cpp
                 */
                factorGraph fg;
                createHMT(fg,allNodeIndices.at(osCounter).data, allParentIndices.at(osCounter).data,summationIndices,0,vSize,indicesXD,type,probO,transProbSeqN.size(),transProbDivN.at(0).size());

                //createHMT( std::vector<std::vector<double> > transProbDiv, std::vector<std::vector<double> > emProb,std::vector<int> observations,std::vector<int> nodeIndices, std::vector<int> parents,std::vector<std::vector<int> > summationIndices,int maxSum,int vSize,std::vector< std::vector<int> > indicesXD,char type)
                std::vector<int>endNodes=fg.endNodes;
                std::vector<int> startNodes;
                startNodes.push_back(1);
                //
                std::vector< std::vector<std::vector <double> > > transProbDivUD;
                transProbDivUD.push_back(transProbDivN);
                transProbDivUD.push_back(arrangeItems(transProbDivN,indicesXD));
                transitionMatrices tM;
                tM.probDiv=transProbDivUD;
                std::vector< std::vector<std::vector<double> > > startProbT;
  
                startProbT.push_back(transpose(startProbN));
                tM.probStart=startProbT;
                std::vector< std::vector<std::vector<double> > > transProbSeqT;
                transProbSeqT.push_back(transProbSeqN);
                tM.probSeq=transProbSeqT;
                //
                
                factorGraph fgN=sumProductAlgorithm(fg,fg.endNodes,startNodes,tM);
                /*
                 if(~isempty(strainParents))
                 tP=tabulate(strainParents(:,1));
                 tP=tP(:,2);
                 else
                 tP=1;
                 end
                 */
                int numberNodes=fgN.factorNodes.size();
                std::vector<double> tP;
                int maxParentValue=0;
                int endNodeCheck=0;
                if(fgN.strainParents.size()>0){
                    tP=tabulate(getColumn(fgN.strainParents,0));
                    //%%%%%%freaky end nodes
                    for(int e=0;e<endNodes.size();e++){
                        int eN=endNodes.at(e);
                        if(!(fgN.factorNodes.at(osCounter).msg_out_down.size()>0)){
                            //if(~isempty(find(isnan(fnN.at(eN).msg_out_down.at(1).msg))))
                            //endNodeCheck=1;
                            //}
                        }
                    }
                    //emptyCells = cellfun('isempty', EijSeq); //%delete the empty cells (multiple observations)
                    //EijSeqT(emptyCells) = [];
                    //%%%%%%%%%%%%%%%%%%%%%%%%%
                    //~isempty(find(isnan(fnN(2).msg_out_up{1}.msg))) ||
                    
                    Max maxParents;
                    maxParents = max(tP);
                    maxParentValue=maxParents.value;
                }
                if(maxParentValue >=3 ||  endNodeCheck ){
                    skipSequences.at(osCounter)=1;
                    continue;
                }else{
                    
                    int N_v=allObservationSequences.at(osCounter).data.size();
                    //emptyCells = cellfun('isempty', EijSeq); //delete the empty cells (multiple observations)
                    //EijSeq(emptyCells) = [];
                    //%%%%%
                    // emptyCells = cellfun('isempty', EijDiv); //delete the empty cells (multiple observations)
                    //
                    eijSeq eijSeqObs;
                    eijSeqObs.data=calculateEijSeq(fgN,transProbSeqN,probO,obsProbSize);
                    allEijSeq.at(osCounter)=eijSeqObs;
                    //
                    eijDiv eijDivObs;
                    eijDivObs.data = calculateEijDiv(fgN,transProbDivN,probO,stateIndicesSingle,obsProbSize,nHStates,nDStates);
                    allEijDiv.at(osCounter)=eijDivObs;
                    //
                    ypsilon ypsilonObs;
                    ypsilonObs.data=calculateYpsilon(fgN,probO,obsProbSize);
                    allYpsilon.at(osCounter)=ypsilonObs;
                    
                }
                //
                double logScale=0;
                for(int e=0;e<fgN.factorNodes.size();e++){
                    //logScale=logScale+fgN.factorNodes.at(e).msg_out_down_scale.at(0);
                }
                for(int e=0;e<fgN.factorNodes.size();e++){
                    if(fgN.factorNodes.at(e).msg_out_down.size()>0){
                    }
                    //logScale=logScale+fgN.factorNodes.at(e).msg_out_down_scale.at(0);
                }
                for(int e=0;e<fgN.factorNodes.size();e++){
                    if(fgN.factorNodes.at(e).msg_out_up.size()>0){
                    }
                    //logScale=logScale+fgN.factorNodes.at(e).msg_out_down_scale.at(0);
                }
                allProbObs.at(osCounter)=(log(sum(multEl(fgN.factorNodes.at(1).msg_out_down.at(0).msg,fgN.factorNodes.at(2).msg_out_up.at(0).msg)))+fgN.factorNodes.at(1).acc_logS_msg_out_down.at(0).msg+fgN.factorNodes.at(2).acc_logS_msg_out_up.at(0).msg);
                // %allProbObs{os}=log(sum(fnO(1).msg_out_down{1}.*fnO(2).msg_out_up{1}))+fnO(1).acc_logS_msg_out_down(1)+fnO(2).acc_logS_msg_out_up(1);
            }
        }
        
        
        if(counterIt>0){
            ll_old=ll;
        }
        ll=0;
        for(int osCounter=0;osCounter<nObs;osCounter++){
            if(allObservationSequences.at(osCounter).data.size()>0 && skipSequences.at(osCounter)==0){
                ll=ll+allProbObs.at(osCounter);
            }
        }
        allLL.push_back(ll);
        allLLSingle.push_back(allProbObs);
        if(counterIt>0){
            if((abs(ll_old-ll)/(abs(ll_old) + abs(ll) + eps))<1e-04 || counterIt>30){
                converged=1;
                break;
            }
        }
        std::vector <double> piN=calculatePI(nHStates,  allYpsilon, skipSequences);
        startProbN.at(0)=piN;
        transProbSeqN=calculateTransProbSeq(nHStates,  skipSequences, allEijSeq );
        //
        transProbDivN=calculateTransProbDiv(nHStates,nDStates, skipSequences,allEijDiv);
        emProbN=calculateBJ(allObservationSequences,nHStates, allYpsilon, skipSequences,oStates);
        /*
         //%%Check first if old model matches convergence criterion%%%%%%%%%%%%%%%
         allLL.push_back(ll);
         ll_old=ll;
         ll=0;
         //ap=zeros(nObs);
         for(int os=0; os <nObs; os++){
         if(sllObservationSequences.at(osCounter).size()>0 && skipSequences.at(os)==0){
         ll=ll+allProbObs.at(os);
         //ap.at(os)=allProbObs.at(os);
         
         }
         
         }
         //allLLSingle.push_back(ap);
         
         if((abs(ll_old-ll)/(abs(ll_old) + abs(ll) + eps))<1e-04 || counterIt>30){
         converged=1;
         break
         }
         */
        cout<<"Log-Likelihood: "<<ll<<endl;
        iterations=iterations+1;
        counterIt=counterIt+1;
    }
    outputD bmOutput;
    bmOutput.transProbSeq=transProbSeqN;
    bmOutput.transProbDiv=transProbDivN;
    bmOutput.prior=startProbN;
    bmOutput.emProb=emProbN;
    bmOutput.allLL=allLL;
    bmOutput.allLLSingle=allLLSingle;
    //printVector(allLL);
    return bmOutput;
}

