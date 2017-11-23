function [ transProbSeqN,transProbDiv,prior,allLL,emProb,mu,sigma ] = HFMfit(features,numStates,indices,priorInit,transmatInit,transInitDiv,type,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% HFMfit Fits a Hidden Factor Graph Model to data.
%
%[ transmat,transmatDiv,prior,allLL,emProb] = HFMfit(features,numStates,indices,priorInit,transmatInit,transInitDiv,type,varargin )
%[ transmat,transmatDiv,prior,allLL,mu,Sigma] = HFMfit(features,numStates,indices,priorInit,transmatInit,transInitDiv,type,varargin )
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
% priorInit(i) = Pr(Q(1) = i), 
% transmatInit -  the initial guess of the transition probabilities of sequential observations, transmatInit(i,j) = Pr(Q(t+1)=j | Q(t)=i)
% transInitDiv -  the initial guess of the transition probabilities of dividing cells, transmatInit(i,j) = Pr(Q(t+1)=j | Q(t)=i)
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
% If the number of mixture components differs depending on Q, just set  the trailing
% entries of mixmat to 0, e.g., 2 components if Q=1, 3 components if Q=2,
% then set mixmat(1,3)=0. In this case, B2(1,3,:)=1.0.
%
%OUTPUT:
%
% maxPath - a N x 4 matrix: Maximum hidden state path (1st column), node
% index (2nd column), branch (of the tree) (2nd column), parent branch (3nd
% column)
%
%EXAMPLE:
%%%%%Fit a discrete treeHFM%%%%%
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
%%%%%Fit a continous treeHFM%%%%%
%%%%%create observation sequences%%%%%
% nHStates = 2;
% O = 2;
% T = 10;
% obs1 = randn(O,10);
% obs2 = randn(O,8);
% data={};
% data{1}=obs1;
% data{2}=obs2;
% %  %%%%% create guesses for prior and transition matrices%%%%%
%  priorInit = ones(nHStates,1)/nHStates;
%  transmatInit = ones(nHStates,nHStates)*(1/nHStates);
%  transInitDiv= ones(nHStates,nHStates*nHStates)*(1/(nHStates*nHStates));
% GMModel = gmdistribution.fit([obs1';obs2'],2);
% muInit=GMModel.mu;
% SigmaInit=GMModel.Sigma;
% %%%%%%create tree topology%%%%%
%  nodeIndices1=[[0,1,2,3,3,4,5,6,7,8]',[1,2,3,4,5,6,7,8,9,10]'];
%  nodeIndices2=[[0,1,2,3,4,4,5,6]',[1,2,3,4,5,6,7,8]'];
%  dataNodeIndices={};
%  dataNodeIndices{1}=nodeIndices1;
%  dataNodeIndices{2}=nodeIndices2;
% %% %%%%% Fit a HFM to the data
% [transProbSeqN,transProbDiv,prior,allLL,emProb,sigma,mu]=HFMfit(data,nHStates,dataNodeIndices,priorInit,transmatInit,transInitDiv,'c','Sigma',SigmaInit,'mu',muInit);

sigmaInit=[];
muInit=[];
emProbN=[];
if(type=='c')
    [sigmaInit,muInit] = process_options(varargin, 'Sigma', [], 'mu', []);
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
parents=cell(1,numel(indices));
nodeIndices=cell(1,numel(indices));
nodeIndicesAll=cell(1,numel(indices));

for entry =1:numel(indices)
    if(~isempty(indices{entry}))
            parents{entry}=indices{entry}(:,2);
            nodeIndices{entry}=indices{entry}(:,1);
            nodeIndicesAll{entry}=[indices{entry}(:,1) indices{entry}(:,2)];
    end
end
%%%%%%delete empty features
emptyCells = cellfun(@isempty,features);
features(emptyCells)=[];
parents(emptyCells)=[];
nodeIndices(emptyCells)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % [ startProbN,transProbSeqN,transProbDivN,emProbN,ll,allLL,muN, SigmaN]= hmt_baumWelchCDMO_CPP_opt( features,priorInit',transmatInit,transInitDiv,[],muInit,sigmaInit,nodeIndicesTemp,summationIndices,stateIndicesSingle,indicesXD,1,sigmaInit)
if(type=='c')
    [transProbSeqN,transProbDiv,prior,allLL,allLLSingle,mu,sigma]= MatlabWrapperBM(priorInit',transmatInit,transInitDiv,emProbN,features,nodeIndices,parents,summationIndices,indicesXD,type,stateIndicesSingle,sigmaInit,muInit);    %[ startProbN2,transProbSeqN2,transProbDivN2,emProbN2,ll2,allLL2,muN2, SigmaN2]= hmt_baumWelchCDMO( features_all,priorInit,transmatInit,transInitDiv,[],muInit,sigmaInit,indicesAll_temp,summationIndices,stateIndicesSingle,0,0,indicesXD )
    %nodeIndicesAll2{1}=[ones(nrows(nodeIndicesAll{:}),1),nodeIndicesAll{:}];
    %[ priorM,transProbSeqNM,transProbDivM,emProbM,ll,allLLM,muM,sigmaM,allEijSeqM,allLLSingleM] =  hmt_baumWelchCDMO_CPP_optLLTest( features,priorInit,transmatInit,transInitDiv,emProbN,muInit,sigmaInit,nodeIndicesAll2,summationIndices,stateIndicesSingle,indicesXD,1,[]);
    %[ startProbMO,transProbSeqMO,transProbDivMO,emProbMO,ll,allLLOld,muMO,SigmaMO]=hmt_baumWelchCDMO_CPP( features,priorInit,transmatInit,transInitDiv,emProbN,muInit,sigmaInit,nodeIndicesAll2,summationIndices,stateIndicesSingle,indicesXD);
    emProb=[];
    allLL=allLL(2:end);
else
    [transProbSeqN,transProbDiv,prior,allLL,emProb]= MatlabWrapperBM(priorInit',transmatInit,transInitDiv,emProbN,parents,features,parents,nodeIndices,summationIndices,indicesXD,type,stateIndicesSingle,sigmaInit,muInit);    %[ startProbN2,transProbSeqN2,transProbDivN2,emProbN2,ll2,allLL2,muN2, SigmaN2]= hmt_baumWelchCDMO( features_all,priorInit,transmatInit,transInitDiv,[],muInit,sigmaInit,indicesAll_temp,summationIndices,stateIndicesSingle,0,0,indicesXD )
end
end

