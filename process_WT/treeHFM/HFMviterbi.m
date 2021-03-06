function [ maxPath,strainParents] = HFMviterbi(features,numStates,indices,prior,transmat,transmatDiv,type,varargin  )
%HFMviterbi= Calculates the most probable hidden state path
%
%Notation: (t) = hidden state, Y(t) = observation, M(t) = mixture variable,
%N = number of observations
%
%INPUTS:
% features -  a cell array of length N with entries features{ex}(:,numStates) if the observations are continous or features{ex} if
% the observations are discrete
% numStates - Number of hidden states of the HFM
% indices - a cell array of length N, the entries give the topology of the
% tree in this observations, they are of the form (parent index, child
% index), e.g. (0 1; 1 2; 1 3)
% prior(i) = Pr(Q(1) = i), 
% transmat -  the initial guess of the transition probabilities of sequential observations, transmatInit(i,j) = Pr(Q(t+1)=j | Q(t)=i)
% transmatDiv -  the initial guess of the transition probabilities of dividing cells, transmatInit(i,j) = Pr(Q(t+1)=j | Q(t)=i)
% type - 'c' for continous observations, 'd' for discrete observations
% 
%
% 
% Parameters for continous observations:
%
% 'Sigma' - Covariance matrix, Sigma(:,:,j,k) = Cov[Y(t) | Q(t)=j, M(t)=k]
% 'mu' -  means of the Gaussian distributions, mu(:,j,k) = E[Y(t) | Q(t)=j, M(t)=k ]
% 'max_iter' - max number of EM iterations [10]
% 'thresh' - convergence threshold [1e-4]
% 'verbose' - if 1, print out loglik at every iteration [1]
% 'cov_type' - 'full', 'diag' or 'spherical' ['full']
%
% To clamp some of the parameters, so learning does not change them:
% 'adj_prior' - if 0, do not change prior [1]
% 'adj_trans' - if 0, do not change transmat [1]
% 'adj_mix' - if 0, do not change mixmat [1]
% 'adj_mu' - if 0, do not change mu [1]
% 'adj_Sigma' - if 0, do not change Sigma [1]
%
%OUTPUT:
%
% maxPath - a N x 4 matrix: Maximum hidden state path (1st column), node
% index (2nd column), branch (of the tree) (2nd column), parent branch (3nd
% column)
%
%EXAMPLE:
%%%%%create observation sequences%%%%%
%  nOStates = 3;
%  nHStates = 2;
%  %%
%  obs1=randi(nOStates,1,10);
%  obs2=randi(nOStates,1,8);
%  data={};
%  data{1}=obs1;
%  data{2}=obs2;
%  %%%%%create tree topology%%%%%
%  nodeIndices1=[[1,2,3,4,5,6,7,8,9,10]',[0,1,2,3,3,4,5,6,7,8]'];
%  nodeIndices2=[[1,2,3,4,5,6,7,8]',[0,1,2,3,4,4,5,6]'];
%  dataNodeIndices={};
%  dataNodeIndices{1}=nodeIndices1;
%  dataNodeIndices{2}=nodeIndices2;
%  %%%%% create guesses for prior and transition matrices%%%%%
%  priorInit = ones(nHStates,1)/nHStates;
%  transmatInit = ones(nHStates,nHStates)*(1/nHStates);
%  transInitDiv= ones(nHStates,nHStates*nHStates)*(1/(nHStates*nHStates));
%  obsmatInit = rand(nHStates,nOStates);
%  % %%%%% Fit a HFM to the data
%  [transProbSeqN,transProbDiv,prior,emProb,allLL]=HFMfit(data,nHStates,dataNodeIndices,priorInit,transmatInit,transInitDiv,'d','observationProb',obsmatInit);
% % %%%%% Calculate the most probable hidden state path
%   [ maxPath] = HFMviterbi(obs1,nHStates,nodeIndices1,prior,transProbSeqN,transProbDiv,'d','observationProb',emProb)
% %%%%% Plot the path
%  [ submatrix ] = drawViterbiTree(maxPath);
%   colors=[1,1,1;1 0 0; 0 1 0; 0 0 1];
%   figure;
%   image(submatrix+1);
%   colormap(colors)




sigma=[];
mu=[];
emProbN=[];
if(type=='c')
    [sigma,mu] = process_options(varargin, 'Sigma', [], 'mu', []);
else
    [emProbN] = process_options(varargin, 'observationProb', []);
end


summationIndices=[];
endState=numStates;
startState=1;
for ns=1:(numStates);
    summationIndices=[summationIndices;startState:endState];
    startState=endState+1;
    endState=endState+numStates;
end
%create state indices
stateIndicesSingle=zeros(numStates*numStates,2);
countS=1;
for nS1=1:numStates
    for nS2=1:numStates
        stateIndicesSingle(countS,:)=[nS1 nS2];
        countS=countS+1;
    end
end
%%%
nHStates=numStates;
indicesXU=1:(numStates*numStates);
indicesXD=[];
for iX1= 1: nHStates
    iX1N=iX1;
    for iX2= 1: nHStates
        indicesXD=[indicesXD,iX1N];
        iX1N=iX1N+nHStates;
    end
end
vSize=6;
%%%%%%%%%
  %%%%%%%%%  //startProb transProbSeqN transProbDiv emProbN nodeIndices parents summationIndices indicesXD stateIndicesSingle type lenObs mu sigma observationSequenceC observationSequenceD dimensionData

%%%%%%%%%
if(strcmp(type,'c'))
    lenObs=size(features,2);
    dimData=size(features,1);
     %startProb transProbSeqN transProbDiv emProbN nodeIndices parents summationIndices indicesXD stateIndicesSingle type lenObs mu sigma observationSequenceC observationSequenceD dimensionData
    [maxPath,strainParents]=MatlabWrapperMS(prior,transmat,transmatDiv,emProbN,indices(:,1)',indices(:,2)',summationIndices,indicesXD,stateIndicesSingle,type,lenObs,sigma,mu,features,dimData);
else
    lenObs=numel(features);
    dimData=size(features,1);
    [maxPath,strainParents]=MatlabWrapperMS(prior,transmat,transmatDiv,emProbN,indices(:,1)',indices(:,2)',summationIndices,indicesXD,stateIndicesSingle,type,mu,sigma,lenObs,features,dimData);
end
%%%%%%%%%

end

