

#ifndef __FGBMA__FactorGraph__
#define __FGBMA__FactorGraph__

#include <stdio.h>
#include "FactorGraph.h"
#include <vector>

struct msgD{
    std::vector<double> msg;
    int index;
};
struct msgB{
    std::vector<int> bd;
    int index;
};
struct msgAccLog{
    double msg;
    int index;
};
struct factorNode {         /* deklariert den Strukturtyp person */
    //
    std::vector<int> f_parents;
    std::vector<int> f_children;
    std::vector<int> v_parents;
    std::vector<int> v_children;
    std::vector<int> v_neighbors;
    //
    int order;
    int leaf;
    int division;
    int nHStates;
    int stateRep;
    int N_v_neighbors;
    int N_msg_in_sofar;
    int N_msg_out_remain;
    int N_msg_in_down;
    int N_msg_in_up;
    int observedFactor;
    //
    std::vector<std::vector<double> > startProb;
    std::vector<std::vector<std::vector<double> > > prob;
    std::vector<std::vector<double> > probO;
    //
    std::vector<int>  indicesXD;
    std::vector<std::vector<int> >  summationIndices;
    std::vector<msgD > msg_out_up;
    std::vector<msgD> msg_out_down;
    std::vector<msgD > msg_in_down;
    std::vector<msgD > msg_in_up;
    //
    std::vector<double> msg_out_up_scale;
    std::vector<double> msg_out_down_scale;
    std::vector<double> msg_in_down_scale;
    std::vector<double> msg_in_up_scale;
    //
    std::vector<msgAccLog> acc_logS_msg_in_up;
    std::vector<msgAccLog> acc_logS_msg_in_down;
    std::vector<msgAccLog> acc_logS_msg_out_up;
    std::vector<msgAccLog> acc_logS_msg_out_down;
    //
    std::vector<msgB> backwardDown;
    std::vector<msgB> backwardUp;
};
struct variableNode {         /* deklariert den Strukturtyp person */
    //
    std::vector<int> f_parents;
    std::vector<int> f_children;
    std::vector<int> v_parents;
    std::vector<int> v_children;
    std::vector<int> f_neighbors;
    std::vector<int> v_neighbors;
    //
    int observed;
    int strain;
    int order;
    int leaf_node;
    int size;
    std::vector<std::vector<double> > value;
    //
    int N_v_neighbors;
    int N_f_neighbors;
    int N_msg_in_sofar;
    int N_msg_out_remain;
    int N_msg_up_remain;
    int  N_msg_down_remain;
    int N_msg_in_down;
    int N_msg_in_up;
    
};
struct factorGraph {
    std::vector<factorNode > factorNodes;
    std::vector<variableNode > variableNodes;
    int maxStrain;
    int maxIndex;
    std::vector<int> singleTransitions;
    std::vector<int> multipleTransitions;
    std::vector<int> endNodes;
    std::vector<int> startNodes;
    std::vector< std::vector<int> > strainParents;
};
struct transitionMatrices{
    std::vector<std::vector<std::vector<double> > > probDiv;
    std::vector<std::vector<std::vector<double> > > probSeq;
    std::vector< std::vector<std::vector<double> > >  probStart;
};

int findMessage(std::vector<msgAccLog> msgV,int index);
int findMessage(std::vector<msgD> msgV,int index);
int findMessage(std::vector<msgB> msgV,int index);
std::vector<int> unionV(std::vector<int>A, std::vector<int>B);
void createHMT(factorGraph &fg, std::vector<int> nodeIndices, std::vector<int> parents,std::vector<std::vector<int> > summationIndices,int maxSum,int vSize,std::vector<int>  indicesXD,char type,std::vector<std::vector<double> > observationProbabilities,int sizeMsgSeq,int sizeMsgDiv);
void addFactorNode(factorNode &f,int index, std::vector<std::vector<double> > probO ,std::vector<int> v_children,std::vector<int> v_parents,std::vector<int> f_children,std::vector<int> f_parents,std::vector<std::vector<int> > summationIndices,int maxNumMessages,int maxSum,int division,int observedFactor,std::vector<int>indicesXD,int size);
void addVariableNode(variableNode &v, int index,int observed,int leaf_node,std::vector<int> v_children,std::vector<int> v_parents,std::vector<int> f_children,std::vector<int> f_parents,int size, int strain);
void drawFactorNode(factorNode f);
void printFactorNode(factorNode f);
void printFactorGraph(struct factorGraph fg);
#endif /* defined(__FGBMA__FactorGraph__) */
