//
//  FactorGraph.cpp
//  FGBMA
//
//  Created by Henrik Failmezger on 16.01.15.
//  Copyright (c) 2015 Henrik Failmezger. All rights reserved.
//

#include "FactorGraph.h"

//
//  SumProduct.cpp
//  FGBMA
//
//  Created by Henrik Failmezger on 16.01.15.
//  Copyright (c) 2015 Henrik Failmezger. All rights reserved.
//


#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include "helpFunctions.h"
#include "time.h"
#include "MultivariateGaussian.h"

using namespace std;





std::vector<int> unionV(std::vector<int>A, std::vector<int>B)
{
    sort(A.begin(),A.end());
    sort(B.begin(),B.end());
    vector<int> v(A.size()+B.size());
    
    vector<int>::iterator it= set_union(A.begin(),A.end(),B.begin(),B.end(),v.begin());
    v.resize(it-v.begin());
    return v;
}





void addFactorNode(factorNode &f,int index, std::vector<std::vector<double> > probO ,std::vector<int> v_children,std::vector<int> v_parents,std::vector<int> f_children,std::vector<int> f_parents,std::vector<std::vector<int> > summationIndices,int maxNumMessages,int maxSum,int division,int observedFactor,std::vector<int>indicesXD,int size){
    f.order=index;
    if(index==1){
        f.leaf=1;
    }else{
        f.leaf=0;
    }
    f.indicesXD=indicesXD;
    f.division=division;
    f.f_parents=f_parents;
    f.f_children=f_children;
    f.v_parents=v_parents;
    f.v_children=v_children;
    
    if(f.v_parents.size()==0 && f.v_children.size()==0){
        std::vector<int> v_neighbors;
        f.v_neighbors=v_neighbors;
    }else if(f.v_parents.size()==0){
        f.v_neighbors=f.v_children;
    }else if(f.v_children.size()==0){
        f.v_neighbors=f.v_parents;
    }else{
        f.v_neighbors=unionV(f.v_parents,f.v_children);
    }
    
    f.N_v_neighbors=f.v_neighbors.size();
    f.N_msg_in_sofar = 0;
    f.N_msg_out_remain = f.N_v_neighbors;
    f.N_msg_in_down =0;
    f.N_msg_in_up = 0;
    f.summationIndices=summationIndices;
    
    //
    //f.msg_out_up_scale=zeros(maxNumMessages);
    //f.msg_out_down_scale=zeros(maxNumMessages);
    //f.msg_in_down_scale=zeros(maxNumMessages);
    //f.msg_in_up_scale=zeros(maxNumMessages);
    //
    
    //f.prob=probS;
    f.probO=probO;
    if(index == 1){
        if(!maxSum){
            msgD msg_in_down;
            msg_in_down.index=1;
            msg_in_down.msg=ones(size);
            f.msg_in_down.push_back(msg_in_down);
            msgAccLog mAL;
            mAL.index=1;
            mAL.msg=log(1.0/sum(ones(size)));
            f.acc_logS_msg_in_down.push_back(mAL);
        }else{
            msgD msg_in_down;
            msg_in_down.index=1;
            msg_in_down.msg=zeros(size);
            f.msg_in_down.push_back(msg_in_down);
            msgAccLog mAL;
            mAL.index=1;
            mAL.msg=log(1.0/sum(zeros(size)));
            f.acc_logS_msg_in_down.push_back(mAL);
        }
    }
    if(index ==maxNumMessages){
        msgAccLog mAL;
        mAL.index=1;
        mAL.msg=log(1.0/sum(ones(size)));
        f.acc_logS_msg_in_down.push_back(mAL);
    }
}
void addVariableNode(variableNode &v, int index,int observed,int leaf_node,std::vector<int> v_children,std::vector<int> v_parents,std::vector<int> f_children,std::vector<int> f_parents,int size, int strain){
    v.v_children=v_children;
    v.v_parents=v_parents;
    v.f_children=f_children;
    v.f_parents=f_parents;
    v.order=index;
    //v.v_neighbors=union(v.v_parents,v.v_children);
    //
    if(v.v_parents.size()==0 && v.v_children.size()==0){
        std::vector<int> v_neighbors;
        v.v_neighbors=v_neighbors;
    }else if(v.v_parents.size()==0){
        v.v_neighbors=v.v_children;
    }else if(v.v_children.size()==0){
        v.v_neighbors=v.v_parents;
    }else{
        v.v_neighbors=unionV(v.v_parents,v.v_children);
    }
    //
    v.N_v_neighbors=v.v_neighbors.size();
    //
    if(v.f_parents.size()==0 && v.f_children.size()==0){
	    	  std::vector<int> f_neighbors;
	    	  v.f_neighbors=f_neighbors;
	   }else if(v.f_parents.size()==0){
           v.f_neighbors=v.f_children;
       }else if(v.f_children.size()==0){
           v.f_neighbors=v.f_parents;
       }else{
           
           v.f_neighbors=unionV(v.f_parents,v.f_children);
           
       }
    
    //
    v.N_f_neighbors = v.f_neighbors.size();
    v.N_msg_in_sofar = 0;
    v.N_msg_out_remain = v.N_f_neighbors;
    v.N_msg_down_remain = v.f_children.size();
    v.N_msg_up_remain = 0;
    v.leaf_node=leaf_node;
    v.size=size;
    v.strain=strain;
    //v(index).size=size_node_list(index);
    
    if(observed){ // j is observed?
        v.observed = 1;
        std::vector <double> valueTMP;
        valueTMP.push_back(1.0);
        valueTMP.push_back(0.0);
        v.value.push_back(valueTMP);
    }else{ //  j is hidden
        v.observed = 0;
    }
}
/*
void drawFactorNode(factorNode f){
    cout << "Index:" <<f.order<<endl;
    cout << "Factor Parents:" <<endl;
    printVector(f.f_parents);
    cout <<endl<< "Factor Childs:" <<endl;
    printVector(f.f_children);
    cout <<endl<< "Variable Parents:" <<endl;
    printVector(f.v_parents);
    cout <<endl<< "Variable Childs:" <<endl;
    printVector(f.v_children);
    
}
void printFactorNode(factorNode f){
    cout << "Index:" <<f.order<<endl;
    cout << "Factor Parents:" <<endl;
    printVector(f.f_parents);
    cout <<endl<< "Factor Childs:" <<endl;
    printVector(f.f_children);
    cout <<endl<< "Variable Parents:" <<endl;
    printVector(f.v_parents);
    cout <<endl<< "Variable Childs:" <<endl;
    printVector(f.v_children);
    
}
void printFactorGraph(struct factorGraph fg){
    cout<<"----------- FACTOR NODES ------------"<<endl;
    for(int i=300; i<400;i++){
        cout<<"Index:"<<i<<endl;
        cout<<fg.factorNodes.at(i).order<<endl;
        cout<<"factor parents:"<<endl;
        printVector(fg.factorNodes.at(i).f_parents);
        cout<<"variable parents:"<<endl;
        printVector(fg.factorNodes.at(i).v_parents);
        cout<<"factor children:"<<endl;
        printVector(fg.factorNodes.at(i).f_children);
        cout<<"variable children:"<<endl;
        printVector(fg.factorNodes.at(i).v_children);
        cout<<"variable neighbors:"<<endl;
        printVector(fg.factorNodes.at(i).v_neighbors);
        cout<<"---"<<endl;
    }
 
}
*/
void createHMT(factorGraph &fg,std::vector<int> nodeIndices, std::vector<int> parents,std::vector<std::vector<int> > summationIndices,int maxSum,int vSize,std::vector<int>  indicesXD,char type,std::vector<std::vector<double> > observationProbabilities,int sizeMsgSeq,int sizeMsgDiv){
   fg.factorNodes.resize(max(nodeIndices).value*2);
    fg.variableNodes.resize(max(nodeIndices).value*2);
    //std::vector<factorNode > f_list(max(nodeIndices).value*2); //factor node list
    //std::vector<variableNode > v_list(max(nodeIndices).value*2); //variable node list
    int lenObs=nodeIndices.size();
    std::vector<std::vector<int> > connectivity;
    for (int o=0; o<lenObs;o++){
               std::vector<int> nodeParent;
               nodeParent.push_back(nodeIndices[o]);
               nodeParent.push_back(parents[o]);
               connectivity.push_back(nodeParent);
    }
    //
    //printVector(connectivity);
    //
    std::vector<std::vector<double> > probO;
    int N=lenObs;
    int maxStrain=0;
    int strain=0;
    std::vector< std::vector<int> > strainParents;
    std::vector<int> endNodes;
    std::vector<int> singleTransitions;
    std::vector<int> multipleTransitions;
    Max maxV=max(nodeIndices);
    int maxIndex=maxV.value;
    int jNode;
    for (int j=0 ; j<lenObs;j++){
        jNode=nodeIndices.at(j);
        //cout<<"jNode: "<<jNode<<endl;
        std::vector<int> f_parents;
        std::vector<int> f_children;
        std::vector<int> v_parents;
        std::vector<int> v_children;
        //
        std::vector<int> f_parents_s;
        std::vector<int> f_children_s;
        std::vector<int> v_parents_s;
        std::vector<int> v_children_s;
        if(j==0){
            //printVector(connectivity);
            //cout<<"found: "<<findIndex(connectivity,1,jNode)[0]<<endl;
            
            std::vector<int> vT=connectivity.at(findIndex(connectivity,1,jNode)[0]); //can only be there 1 time?
            if(jNode==1){
            }
            v_children_s.push_back(jNode);
            v_children.push_back(vT.at(0));
            f_parents.push_back(1);
            //
            std::vector<int> f_childrenT;
            f_childrenT.push_back(vT[0]);
            f_childrenT.push_back(jNode+maxIndex);
            //
            //cout<<"prob start"<<endl;
            std::vector<std::vector<double> > probOStart(max(nodeIndices).value*2);;
            /*
            if(type=='d'){
                probOStart.at(j)=emProb.at(observations.at(j-1));
            }else{
                probOStart.at(j)=observationProbabilities.at(j-1);  //////BEI GELEGENHEIT NOCHMAL CHECKEN!!!!!
            }
            */
            
            probOStart.at(jNode)=observationProbabilities.at(jNode-1);
            //cout<<"probO"<<endl;
            //v_children.push_back(j+maxIndex);
            //int index,std::vector<double> probS,std::vector<double> probO ,std::vector<int> v_children,std::vector<int> v_parents,std::vector<int> f_children,std::vector<int> f_parents,std::vector<std::vector<int> > summationIndices,int maxNumMessages,int maxSum,int division,int observedFactor,std::vector<int>indicesXD
            int maxNodeIndex=max(nodeIndices).value*2;
            factorNode f;
            addFactorNode(f,jNode,probOStart,v_children_s,v_parents_s,f_childrenT,f_parents_s,summationIndices,maxNodeIndex,maxSum,0,0,indicesXD,sizeMsgSeq);
            //int observedFactor,std::vector<int>indicesXD)
            //f_list.at(j)=f;
            fg.factorNodes.at(jNode)=f;
            //variable node
        }else{
            std::vector<int> indicesChildren=findIndex(connectivity,1,jNode);
            for(int i=0; i<indicesChildren.size();i++){
                v_children.push_back(connectivity.at(indicesChildren.at(i))[0]);
            }
            
            std::vector<int> indicesParents=findIndex(connectivity,0,jNode);
            if(indicesParents.size() >0){ // if a parent was found
                v_parents.push_back(connectivity.at(indicesParents[0]).at(1)); //only 1 allowed
                f_parents.push_back(v_parents.at(0)+1);
            }
            //check for common parent
            if(v_parents.size()>0){
                if(countIndex(connectivity,1,v_parents.at(0))>1){
                    strain=maxStrain+1; //add new strain
                    maxStrain=strain;
                    std::vector<int> strainParent;
                    
                    strainParent.push_back(fg.variableNodes.at(v_parents.at(0)).strain);
                    strainParent.push_back(strain);
                    strainParents.push_back(strainParent);
                }else{
                    strain=fg.variableNodes.at(v_parents.at(0)).strain  ;
                }
            }else{
                strain=maxStrain;
            }
        }
        if(v_children.size() >0){
            f_children.push_back(jNode+1);
        }else{
            endNodes.push_back(jNode);
            std::vector<int> f_children;
        }
        int division=0;
        int MsgSizeTemp=sizeMsgSeq;
        if(v_children.size()>1){
            division=1;
            int MsgSizeTemp=sizeMsgDiv;
        }
        variableNode v;
        addVariableNode(v,jNode,0,0,v_children,v_parents,f_children,f_parents,vSize,strain);
        //v_list.at(j)=v;
        fg.variableNodes.at(jNode)=v;
        int leaf_node=0;
        if(v_children.size()>0){
            std::vector<int> factorNode_f_children=add(v_children,1);
            probO.clear();
            //cout<<"probO"<<endl;
            std::vector<std::vector<double> > probO(max(nodeIndices).value*2);
            /*
            if(type=='d'){
                for(int i=0; i<v_children.size(); i++){
                    int v_child=v_children.at(i);
                    probO.at(v_child)=emProb.at(observations.at(v_child-1));
                }
            }else{
                for(int i=0; i<v_children.size(); i++){
                    int v_child=v_children.at(i);
                    probO.at(v_child)=observationProbabilities.at(v_child-1);
                }
            }
             */
            for(int i=0; i<v_children.size(); i++){
                int v_child=v_children.at(i);
                probO.at(v_child)=observationProbabilities.at(v_child-1);
            }
            std::vector<int> v_parentsF;
            v_parentsF.push_back(jNode);
            //(j,startProb,probO,v_children,v_parents,f_childrenT,f_parents,summationIndices,max(nodeIndices)*2,maxSum,0,0,indicesXD);
            int maxNodeIndex=max(nodeIndices).value*2;
            factorNode f;
            //cout<<"addFactorNode"<<endl;
            addFactorNode(f,jNode+1,probO,v_children,v_parentsF,factorNode_f_children,f_parents,summationIndices,maxNodeIndex,maxSum,division,0,indicesXD,MsgSizeTemp);//factor node right
            if(v_children.size()>1){
                multipleTransitions.push_back(jNode+1);
            }else{
                singleTransitions.push_back(jNode+1);
            }
            fg.factorNodes.at(jNode+1)=f;
            
        }
        
	   }
    
    //fg.factorNodes=f_list;
    //fg.variableNodes=v_list;
    fg.maxStrain=maxStrain;
    fg.maxIndex=maxIndex;
    fg.multipleTransitions=multipleTransitions;
    fg.singleTransitions=singleTransitions;
    fg.endNodes=endNodes;
    fg.strainParents=strainParents;
    //cout<<"factor graph created"<<endl;
}


int findMessage(std::vector<msgAccLog> msgV,int index){
    int indF=-1;
    for(int i=0; i<msgV.size(); i++){
        if(msgV.at(i).index==index){
            return i;
        }
    }
    return indF;
}
int findMessage(std::vector<msgB> msgV,int index){
    int indF=-1;
    for(int i=0; i<msgV.size(); i++){
        if(msgV.at(i).index==index){
            return i;
        }
    }
    return indF;
}
int findMessage(std::vector<msgD> msgV,int index){
    int indF=-1;
    for(int i=0; i<msgV.size(); i++){
        if(msgV.at(i).index==index){
        return i;
        }
    }
    return indF;
}