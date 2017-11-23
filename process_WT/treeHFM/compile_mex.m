mex -g  src/MatlabWrapperBM.cpp src/SumProduct.cpp  src/BaumWelch.cpp src/FactorGraph.cpp src/helpFunctions.cpp  src/MultivariateGaussian.cpp src/ParamContainerEmissions.cpp src/EmissionFunction.cpp src/matUtils.cpp src/MemoryAllocation.cpp -Dchar16_t=uint16_t -lmwlapack

mex  src/MatlabWrapperMS.cpp src/MaxSum.cpp  src/SumProduct.cpp  src/FactorGraph.cpp src/helpFunctions.cpp  src/MultivariateGaussian.cpp src/ParamContainerEmissions.cpp src/EmissionFunction.cpp src/MemoryAllocation.cpp src/matUtils.cpp  -Dchar16_t=uint16_t -lmwlapack
