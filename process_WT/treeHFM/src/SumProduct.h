

#ifndef __FGBMA__SumProduct__
#define __FGBMA__SumProduct__

#include <stdio.h>
#include "helpFunctions.h"
#include "FactorGraph.h"

void sendMessageFToVHMT( factorGraph &fg, int j, int i,transitionMatrices tM);
void sendMessageVToFHMT( factorGraph &fg, int j, int i);
factorGraph sumProductAlgorithm(factorGraph fgS, std::vector<int> endNodes, std::vector<int> startNodes, transitionMatrices tM);




#endif /* defined(__FGBMA__SumProduct__) */
