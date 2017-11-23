

#ifndef __FGBMA__MaxSum__
#define __FGBMA__MaxSum__
#include <queue>
#include "FactorGraph.h"
#include "helpFunctions.h"
void sendMessageVToFMSHMT( factorGraph &fg, int j, int i);
void sendMessageFToVMSHMT( factorGraph &fg, int j, int i,transitionMatrices tM);
std::vector< std::vector<int> > maxSumAlgorithm(factorGraph fgS, std::vector<int> endNodes, std::vector<int> startNodes, transitionMatrices tM, int N_V);
#include <stdio.h>

#endif /* defined(__FGBMA__MaxSum__) */
