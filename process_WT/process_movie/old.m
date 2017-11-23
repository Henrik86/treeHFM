%%%%%%%%%%%%%%%%%%%%%%%%%%
transProbDivND=transProbDivN+transProbDivN(:,indicesXD)
transProbDivND=transProbDivND./repmat(sum(transProbDivND,2),1,ncols(transProbDivND))

%transProbDivDoubled=zeros(nrows(transProbDivN),nrows(indicesSingleToDoubleN(:,2)));
%transProbDivDoubled(:,indicesSingleToDouble(:,2))= transProbDivN(:,indicesSingleToDoubleN(:,1));
%transProbDivUD{1}=transProbDivDoubled;
%transProbDivUD{2}=transProbDivDoubled(:,indicesXD);
tic
startProbM=startProbN;
transProbSeqM=transProbSeqN;
transProbDivM=transProbDivN;%Small factor added to avoid zeros
emProbM=emProbN;
transProbDivUD{1}=transProbDivND;
transProbDivUD{2}=transProbDivND(:,indicesXD);
%transProbDivUD{1}=transProbDivND;
%transProbDivUD{2}=transProbDivND(:,indicesXD);
maxSum=1;
space=3
numFrames1=length(m1Data.handles.Measurements.Image.Count_Nuclei);
numCells1=max(unique(trajIndex1));
stateMatrix1=zeros(numCells1*3+numCells1*space,numFrames1);
treeCounter=1
indexStates=ncols(cellClusterM1)+1;
cellClusterM1(:,indexStates)=zeros(1,nrows(cellClusterM1));
indexStrain=ncols(cellClusterM1)+2;
cellClusterM1(:,indexStrain)=zeros(1,nrows(cellClusterM1));
trajStrainParents={}%assignment strain/parent strain per Traj
trajCellStrain={}
for trajInd1=1:length(trajIndices1)
    %if length(find(ismember(indicesAll{trajIndices1(trajInd1)}(:,1),1)))<=2
    if nrows(indicesAll{trajIndices1(trajInd1)}(:,1))>1
        %[f1,v1]=createHMT(
        %startProbN',transProbSeqN,transProbDivUD,emProbN,centersTraj_all{trajIndices1(trajInd1)},indicesAll{trajIndices1(trajInd1)}(:,2),indicesAll{trajIndices1(trajInd1)}(:,3),summationIndices,maxSum,vSize)
        B = mixgauss_prob(featuresAll{trajIndices1(trajInd1)}, muN, SigmaN);
        [f1,v1,singleTransitions,multipleTransitions,endNodes,strainParents]=createHMTPruning( startProbM,transProbSeqM,transProbDivUD,B',centersTraj_all{trajIndices1(trajInd1)},indicesAll{trajIndices1(trajInd1)}(:,1),indicesAll{trajIndices1(trajInd1)}(:,2),summationIndices,maxSum,vSize,indicesXD,'continous');
        if(~isempty(strainParents))
            t=tabulate(strainParents(:,1))
            t=t(:,2)
        else
            t=1
        end
        if(max(t)<3) %only division in two cells allowed
            N=length(centersTraj_all{trajIndices1(trajInd1)});
            N_v=N;
            N_f=N_v;
            [ f,v,maxPath]= maxSumAlgorithmHMT_Fast(f1,v1,N_v,N_v,endNodes);

            cellStrain=zeros(nrows(maxPath),2);
            cellStrain(maxPath(:,1),1)=maxPath(:,2); %order the entries in the right direction
            cellStrain(maxPath(:,1),2)=maxPath(:,3); %order the entries in the right direction
            trajCellStrain{trajIndices1(trajInd1)}=cellStrain;
            %strainsTrack=unique(cellStrain(:,2))
            %M=zeros(length(strainstrack),length(centersTraj_all{trajIndices1(trajInd1)}));
            M=drawTree(strainParents,'states',cellStrain)
            trajStrainParents{trajIndices1(trajInd1)}=strainParents
            %imshow(M+1,[1,1,1;COLORS_STATES])
            %M=zeros(3,length(centersTraj_all{trajIndices1(trajInd1)}));
            %cellsStrain0=cellStrain(cellStrain(:,2)==0,1);
            %cellsStrain1=cellStrain(cellStrain(:,2)==1,1);
            %cellsStrain2=cellStrain(cellStrain(:,2)==2,1);
            %M(2,1:nrows(cellsStrain0))=cellsStrain0;
            %M(1,(nrows(cellsStrain0)+1):(nrows(cellsStrain0)+1+nrows(cellsStrain1)-1))=cellsStrain1;
            %M(3,(nrows(cellsStrain0)+1):(nrows(cellsStrain0)+1+nrows(cellsStrain2)-1))=cellsStrain2;
            stateMatrix1(treeCounter:(treeCounter+nrows(M)-1),1:ncols(M))=M ; 
            treeCounter=treeCounter+nrows(M)+space;
            if(~isempty(maxPath))
                cellClusterM1(trajIndex1==trajIndices1(trajInd1),indexStates)= cellStrain(:,1);
                cellClusterM1(trajIndex1==trajIndices1(trajInd1),indexStrain)= cellStrain(:,2);
            end
        end
    end
end