function [ submatrix ] = drawViterbiTree(maxPath  )
%INPUTS: The output of HFMviterbi.
%OUTPUT:
%
%The Viterbi path arranged in a matrix.
%
%EXAMPLE:
% nOStates = 3;
% nHStates = 2;
%
% obs=randi(3,1,5); %observation sequence
% data={};
% dataNodeIndices={};
% data{1}=obs;
% nodeIndices=[[1 2 3 4 5]',[0 1 2 2 4]']; %tree structure
% dataNodeIndices{1}=nodeIndices;
% priorInit = ones(nHStates,1)/nHStates;
% transmatInit = ones(nHStates,nHStates)*(1/nHStates);
% transInitDiv= ones(nHStates,nHStates*nHStates)*(1/(nHStates*nHStates))%rr rg rb gr gg gb br bg bb
% obsmatInit = rand(nHStates,nOStates);
% obsmatInit=obsmatInit./repmat(sum(obsmatInit,2),1,nOStates);
% %%%%% Fit a HFM to the data
% [transProbSeqN,transProbDiv,prior,emProb,allLL]=HFMfit(data,nHStates,dataNodeIndices,priorInit,transmatInit,transInitDiv,'d','observationProb',obsmatInit);
% %%%%% Calculate the most probable hidden state path
%  [ maxPath] = HFMviterbi(obs,nHStates,nodeIndices,prior,transProbSeqN,transProbDiv,'d','observationProb',emProb)
% %%%%% Plot the path
% [ submatrix ] = drawViterbiTree(maxPath);
%  colors=[1,1,1;1 0 0; 0 1 0; 0 0 1];
%  figure;
%  image(submatrix+1);
%  colormap(colors)



cellsStrain=zeros(nrows(maxPath),2);
cellsStrain(maxPath(:,1),1)=maxPath(:,2); %order the entries in the right direction
cellsStrain(maxPath(:,1),2)=maxPath(:,3);
%%%%%%%%%%%
parents=[];
for st=maxPath(:,3)'
    parentStrain=maxPath(maxPath(:,3)==st,4);
    stP=unique(parentStrain);
    indexFP=find(stP ~=st);
    if(~isempty(indexFP))
        parents=[parents;stP(indexFP),st]
    end
end
%%%
if(isempty(parents))
           submatrix=cellsStrain(:,1)';
else
    [C,IA,IC] = unique(parents(:,2));
    parents=parents(IA,:)
    nodes=[unique(parents(:))'];
    visitedNodes=zeros(length(nodes),1);
    startNodes=[1];
    matrix=[];
    stopDrawing=0;
    nodeSubtree=cell(1,length(nodes));
    while ~stopDrawing
        for subtree=nodes
            children=parents(parents(:,1)==subtree,2);
            if isempty(children)
                matrix=cellsStrain(cellsStrain(:,2)==subtree,1)';
               nodeSubtree{subtree+1}=matrix;
               visitedNodes(subtree+1)=1;
            else
                child1=children(1);
                child2=children(2);
                if(visitedNodes(child1+1)==1 && visitedNodes(child2+1)==1)
                    %join
                    rowsChild1=nrows(nodeSubtree{child1+1});
                    rowsChild2=nrows(nodeSubtree{child2+1});
                    colsChild1=ncols(nodeSubtree{child1+1});
                    colsChild2=ncols(nodeSubtree{child2+1});
                    matrixNew=cellsStrain(cellsStrain(:,2)==subtree,1)';
                    matrixJoin=zeros(rowsChild1+rowsChild2+nrows(matrixNew),max(colsChild1,colsChild2));
                    matrix2=zeros(rowsChild1+rowsChild2+nrows(matrixNew),ncols(matrixNew));
                    matrixJoin(1:rowsChild1,1:colsChild1,:)=nodeSubtree{child1+1};
                    matrixJoin((rowsChild1+nrows(matrixNew)+1):(rowsChild1+1+nrows(matrixNew)+(rowsChild2-1)),1:colsChild2,:)=nodeSubtree{child2+1};
                    %add the new strain
                    matrix2((rowsChild1+1):(rowsChild1+1+nrows(matrixNew)-1),1:ncols(matrixNew),:)=matrixNew;
                    %
                    nodeSubtree{subtree+1}=[matrix2,matrixJoin];
                    visitedNodes(subtree+1)=1;
                end
            end
            if(isempty(find(visitedNodes==0)))
                stopDrawing=1;
                break
            end
        end
    end
    submatrix=nodeSubtree{1};
end
end

