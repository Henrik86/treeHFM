function [cellClusterM1,featuresAll,indicesAll,trajIndex1,cellCoordinates ] = dataPreprocessing(processedData1,measNames1,featureNames)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
measNames=measNames1;
cellClusterM1=processedData1;


indexGroup1=cellClusterM1(:,ismember(measNames, 'Image_number'));
x_loc1= cellClusterM1(:,ismember(measNames, 'Location_Center_X'));
y_loc1= cellClusterM1(:,ismember(measNames, 'Location_Center_Y'));
cellNumber1=cellClusterM1(:,ismember(measNames, 'Number_Object_Number'));
divisionIndex1=cellClusterM1(:,ismember(measNames, 'DivisionEvent'));
trajIndex1=cellClusterM1(:,ismember(measNames, 'trajIndex'));
cellIndices=[cellClusterM1(:,ismember(measNames,'trajCellNumber')),cellClusterM1(:,ismember(measNames,'trajParentCellNumber' )),divisionIndex1];
%%%
cellCoordinates=[x_loc1,y_loc1,cellNumber1,indexGroup1];

cellClusterM1=cellClusterM1(:,ismember(measNames,featureNames));
[pc,score,latent,tsquare]=princomp(zscore(cellClusterM1));
cm=cumsum(latent)./sum(latent);
indComp=find(cm>=0.90);
cellClusterM1=score(:,1:indComp(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
featuresAll={};
trajIndices1=unique(trajIndex1);
indicesAll={};
trajectories1={};
trajIndices1=unique(trajIndex1);
for trajInd1=1:length(trajIndices1)
    featuresAll{trajIndices1(trajInd1)}=cellClusterM1(trajIndex1==trajIndices1(trajInd1),:)';
    indicesAll{trajIndices1(trajInd1)}=cellIndices(trajIndex1==trajIndices1(trajInd1),:);
end

%%%%%%%%%%%%%%%%%%%PRUNE TRACKS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%assing strain to every cell
trajConnectivity=cell(1,max(trajIndices1));
indicesToExclude=[];
indicesAll_pruned={};
for trajInd1=1:length(trajIndices1)
    trajectory=indicesAll{trajIndices1(trajInd1)};
    connectivity=[];
    for o=1:nrows(trajectory)
       if(o>1)
            connectivity=[connectivity;trajectory(o,1) trajectory(o,2)];
       end
    end
    trajConnectivity{trajIndices1(trajInd1)}=connectivity;
    maxStrain=0;
    strain=0;
    strainParentChild=[];
    for cIndex=1:nrows(trajectory)
        cellIndex=trajectory(cIndex,1);
        if cIndex==1
            trajectory(cIndex,4)=strain;
        else
            v_parents=connectivity(connectivity(:,1)==cellIndex,2)';
            if(length(connectivity(connectivity(:,2)==v_parents))>1)
                neighbor=setdiff(connectivity(connectivity(:,2)==v_parents),j);
                
                strain=maxStrain+1; %add new strain
                maxStrain=strain;
                trajectory(cIndex,4)=strain;
                %
                v_parent=connectivity(connectivity(:,1)==cellIndex,2)';
                v_parentIndex=find(trajectory(:,1)==v_parent);
                v_parentStrain=trajectory(v_parentIndex,4);
                strainParentChild=[strainParentChild;[v_parentStrain,strain]];
                trajectory(v_parentIndex,3)=-1; %add the number at the end of the strain, in order to see if there is a division event afterwars
            else%add strain of parent
                v_parent=connectivity(connectivity(:,1)==cellIndex,2)';
                v_parentIndex=find(trajectory(:,1)==v_parent);
                strainParent=trajectory(v_parentIndex,4);
                trajectory(cIndex,4)=strainParent;
            end
            
        end
        
    end
    %%%%%%%%%%prune strains%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indicesAll{trajIndices1(trajInd1)}=trajectory;
    if(nrows(trajectory)<=1)%Prune very short sequences
        indicesToExclude=[indicesToExclude,trajInd1];
    end
end
%calculate mean cell cycle time
cc_dur_1=[];
allDivisionTimes1=0;
allLifetimes1=0;
for trajInd1=1:length(trajIndices1)
    trajectory=indicesAll{trajIndices1(trajInd1)};
    strains=unique(trajectory(:,4));
    if(length(strains)>=2)
        for strain=strains'
            division=trajectory(find(trajectory(:,4)==strain),1);
            if(length(unique(division))>=3)
                allDivisionTimes1=allDivisionTimes1+1;
                allLifetimes1=allLifetimes1+nrows(division)
                cc_dur_1=[cc_dur_1,ncols(division)];%%length of one strain == length of one division
            end
        end
    end
end
dur1_time=allLifetimes1/allDivisionTimes1;


   %%insert Cell Cycle time in the trajectorie
    cc_dur_1=[];
    cc_time_1={};
    for trajInd1=1:length(trajIndices1)  
        trajectory=indicesAll{trajIndices1(trajInd1)};
        cc_time_traj=zeros(nrows(trajectory),1);
        strains=unique(trajectory(:,4));
        if(length(strains)>1)
            for strain=strains'
                strainIndices=find(trajectory(:,4)==strain);
                nrows(strainIndices)
                cellDivision=trajectory(find(trajectory(:,4)==strain),1);
                lifetime=nrows(cellDivision);
                if(cellDivision(1)==1&&length(find(ismember(cellDivision,-1))) >=1 )
                    timeBetweenDiv=(1:lifetime)/lifetime;
                    cc_time_traj(strainIndices)=timeBetweenDiv;
                elseif(cellDivision(1)==1)%after division to end
                    if(lifetime<=dur1_time)
                        restTime=(1:lifetime)/dur1_time;
                        cc_time_traj(strainIndices)=restTime;
                    else
                        restTime=(1:lifetime)/lifetime;
                        cc_time_traj(strainIndices)=restTime;
                    end
                else%from start to division
                    if(lifetime>dur1_time)
                        timeBeforeDiv=(1:lifetime)/lifetime;
                        cc_time_traj(strainIndices)=timeBeforeDiv;
                    else
                        timeBeforeDiv=((dur1_time-lifetime):(dur1_time-1))/dur1_time;
                        cc_time_traj(strainIndices)=timeBeforeDiv;
                    end
                end
            end
            cc_time_1{trajIndices1(trajInd1)}=cc_time_traj;
            indicesAll{trajIndices1(trajInd1)}=[indicesAll{trajIndices1(trajInd1)},cc_time_traj];
        else
            %no division, label from start to end
            lifetime=nrows(trajectory)
            cc_time_1{trajIndices1(trajInd1)}=(1:lifetime)/lifetime;
            indicesAll{trajIndices1(trajInd1)}=[indicesAll{trajIndices1(trajInd1)},((1:lifetime)/lifetime)']
        end
    end
    
 
    %add the cc timepoints to the data in order to be able to process it -  KD1
    indexCC_time=ncols(cellClusterM1)+1;

    for trajInd1=1:length(trajIndices1)
        trajIndices1(trajInd1);
        %do only use trajectories of a certain duration (else failures)
            size(cellClusterM1(trajIndex1==trajIndices1(trajInd1),:));
            size(cc_time_1{trajIndices1(trajInd1)});
            cellClusterM1(trajIndex1==trajIndices1(trajInd1),indexCC_time)= cc_time_1{trajIndices1(trajInd1)};
    end   
end

