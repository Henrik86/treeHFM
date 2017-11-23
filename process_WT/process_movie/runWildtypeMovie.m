


clear
%path to movie output file
filename1='APC-135395-Neg-N';%
%m1Data=load('/Users/henrikfailmezger/Documents/PHD/Time-Lapse-Movies/Movies/Tracked/LT0042_09--ex2006_03_03--sp2005_05_18--tt17--c5/P1/015--02--03--(3,1)--103860--Neg1-gfp.avi/015--02--03--(3,1)--103860--Neg1-gfp.avi.mat');
m1Data=load('../data/processed/015--02--03--(3,1)--103860--Neg1-gfp.avi.mat');
%
pathToMovie='../data/raw/015--02--03--(3,1)--103860--Neg1-gfp.avi';
%%Create the output directory
outputFolderHMM=strcat('/Users/henrikfailmezger/Documents/PHD/Time-Lapse-Movies/Processed_Movies/Manuskript_ImagesAD/');
mkdir(outputFolderHMM)
%Track the cells over time
[processedData1,measNames1,nucleiCells1]=TrackNucleiOverOneSiteHMT(m1Data,outputFolderHMM,1,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
featureNames={ 'AreaShape_Zernike_6_0'  'AreaShape_Zernike_6_2'  'AreaShape_Zernike_6_4'  'AreaShape_Zernike_6_6'  'Texture_InfoMeas1_Nucleus_3_0'  'AreaShape_Compactness'  'RadialDistribution_RadialCV_Nucleus_2of4'  'AreaShape_Zernike_2_0'  'AreaShape_Zernike_2_2'  'RadialDistribution_FracAtD_Nucleus_2of4'  'AreaShape_Zernike_4_4'  'AreaShape_Zernike_4_2'  'AreaShape_Zernike_4_0'  'AreaShape_Zernike_8_8'  'AreaShape_Zernike_5_1'  'AreaShape_Zernike_8_2'  'AreaShape_Zernike_8_0'  'AreaShape_Solidity'  'AreaShape_Zernike_8_4'  'AreaShape_Eccentricity'  'RadialDistribution_RadialCV_Nucleus_4of4'  'RadialDistribution_FracAtD_Nucleus_4of4'  'AreaShape_MaxFeretDiameter'  'AreaShape_Zernike_7_1'  'AreaShape_Zernike_7_3'  'AreaShape_Zernike_7_5'  'AreaShape_Zernike_7_7'  'AreaShape_EulerNumber'  'Intensity_StdIntensityEdge_Nucleus'  'AreaShape_Zernike_0_0'  'Intensity_MinIntensityEdge_Nucleus'  'AreaShape_MinorAxisLength'  'Texture_Entropy_Nucleus_3_0'  'RadialDistribution_RadialCV_Nucleus_1of4'  'AreaShape_Zernike_3_1'  'AreaShape_Zernike_3_3'  'RadialDistribution_MeanFrac_Nucleus_1of4'  'RadialDistribution_RadialCV_Nucleus_3of4'  'AreaShape_Zernike_5_5'  'AreaShape_Zernike_5_3'  'AreaShape_Perimeter'  'AreaShape_FormFactor'  'RadialDistribution_MeanFrac_Nucleus_2of4'  'AreaShape_Zernike_9_3'  'AreaShape_Zernike_9_1'  'AreaShape_Zernike_9_7'  'AreaShape_Zernike_9_5'  'Intensity_StdIntensity_Nucleus'  'AreaShape_Zernike_9_9'  'Texture_Contrast_Nucleus_3_0'  'RadialDistribution_MeanFrac_Nucleus_3of4'  'Texture_SumEntropy_Nucleus_3_0'  'Intensity_MedianIntensity_Nucleus'  'Texture_DifferenceVariance_Nucleus_3_0'  'AreaShape_Zernike_1_1'  'RadialDistribution_FracAtD_Nucleus_1of4'  'Texture_DifferenceEntropy_Nucleus_3_0'  'RadialDistribution_FracAtD_Nucleus_3of4'  'Intensity_MinIntensity_Nucleus'  'Intensity_IntegratedIntensityEdge_Nucleus'  'Intensity_MaxIntensityEdge_Nucleus'  'Intensity_UpperQuartileIntensity_Nucleus'  'Intensity_MassDisplacement_Nucleus'  'Intensity_MeanIntensity_Nucleus'  'Intensity_MaxIntensity_Nucleus'  'Texture_Correlation_Nucleus_3_0'  'Texture_AngularSecondMoment_Nucleus_3_0'  'Texture_SumVariance_Nucleus_3_0'  'Texture_Variance_Nucleus_3_0'  'AreaShape_Area'  'AreaShape_MinFeretDiameter'  'Intensity_LowerQuartileIntensity_Nucleus'  'Texture_SumAverage_Nucleus_3_0'  'Texture_InfoMeas2_Nucleus_3_0'  'AreaShape_MedianRadius'  'AreaShape_Extent'  'Texture_InverseDifferenceMoment_Nucleus_3_0'  'RadialDistribution_MeanFrac_Nucleus_4of4'  'AreaShape_MajorAxisLength'  'AreaShape_MeanRadius'  'Texture_Gabor_Nucleus_3'  'AreaShape_MaximumRadius'  'AreaShape_Zernike_8_6'  'Intensity_IntegratedIntensity_Nucleus'  'Intensity_MeanIntensityEdge_Nucleus' 
};
[cellClusterM1,featuresAll,indicesAll,trajIndex1,cellCoordinates] = dataPreprocessing(processedData1,measNames1,featureNames);
%%%%%%%%%%%%%%%%%%%%%%%%%%    
COLORS_STATES = [255,0,0;
            0,255,0;
            0,0,255;
            255,255,0;
            255,0,255;
            0,255,255; 
            139,71,38;
            169,169,169;
            160,32,240;
            0,100,0;
            255,193,193;
            255,140,105;
            205,175,149;
            238,118,0;
            205,205,0;
            47,79,79;
            0,0,128;
            56,142,142;
            142,142,56;
            198,113,113;
            92,92,92
            ]/255;%%%COLORS for the Hidden States
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trajIndices1=unique(trajIndex1);
indexCC_time=ncols(cellClusterM1);
numStates=6;
numClusterCenters=numStates;
clusters = kmeans(cellClusterM1(:,1:(indexCC_time-1)),numClusterCenters,'EmptyAction','singleton','replicates',5);
muInit=[];
sigmaInit=[];
counter=1;
for x1 = 1:numClusterCenters
    %indexCC_time has to be the last index
    dataState=cellClusterM1(clusters==x1,1:(indexCC_time-1));
    [mu0, Sigma0] = mixgauss_init(1, dataState', 'full');
    muInit=[muInit,mu0];
    sigmaInit(:,:,counter)=Sigma0;
    counter=counter+1;
end
priorInit = ones(numStates,1)*1/numStates;
transmatInit = ones(numStates,numStates)*(1/numStates);
transInitDiv= ones(numStates,numStates*numStates)*(1/(numStates*numStates));
sigmaInit=sigmaInit+0.08;
[ transProbSeq,transProbDiv,prior,allLL,emProb,mu,sigma ] = HFMfit(featuresAll,numStates,indicesAll,priorInit,transmatInit,transInitDiv,'c','sigma',sigmaInit,'mu',muInit );
%%%% mapping from
indicesXD=[];
for iX1= 1: numStates
    iX1N=iX1;
    for iX2= 1: numStates
        indicesXD=[indicesXD,iX1N];
        iX1N=iX1N+numStates;
    end
end
plotTransitionMatrices(transProbSeq,transProbDiv,COLORS_STATES,indicesXD);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% N=6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create the mapping from first order to second order transitions
stateIndicesSingle=zeros(numStates*numStates,2);
countS=1;
for nS1=1:numStates
    for nS2=1:numStates
        stateIndicesSingle(countS,:)=[nS1 nS2];
        countS=countS+1;
    end
end
mappingMatrix=zeros(numStates,numStates);
for si=1:nrows(stateIndicesSingle)
    cords=stateIndicesSingle(si,:);
    mappingMatrix(cords(1),cords(2))=si;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxSum=1;
space=3;
numFrames1=length(m1Data.handles.Measurements.Image.Count_Nuclei);
numCells1=max(unique(trajIndex1));
stateMatrix1=zeros(numCells1*3+numCells1*space,numFrames1);
treeCounter=1;
indexStates=ncols(cellClusterM1)+1;
cellClusterM1(:,indexStates)=zeros(1,nrows(cellClusterM1));
indexStrain=ncols(cellClusterM1)+2;
cellClusterM1(:,indexStrain)=zeros(1,nrows(cellClusterM1));
trajStrainParents={};
%assignment strain/parent strain per Traj
trajCellStrain={};
for trajInd1=1:length(trajIndices1)
    if nrows(indicesAll{trajIndices1(trajInd1)}(:,1))>1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ maxPath,strainParents] = HFMviterbi(featuresAll{trajIndices1(trajInd1)},numStates,indicesAll{trajIndices1(trajInd1)},prior,transProbSeq,transProbDiv,'c','Sigma',sigma,'mu', mu);
        if(~isempty(maxPath))
            cellStrain=zeros(nrows(maxPath),2);
            cellStrain(maxPath(:,1),1)=maxPath(:,2); %order the entries in the right direction
            cellStrain(maxPath(:,1),2)=maxPath(:,3); %order the entries in the right direction
            M=drawViterbiTree(maxPath);
            stateMatrix1(treeCounter:(treeCounter+nrows(M)-1),1:ncols(M))=M ; 
            treeCounter=treeCounter+nrows(M)+space;
            cellClusterM1(trajIndex1==trajIndices1(trajInd1),indexStates)= cellStrain(:,1);
            cellClusterM1(trajIndex1==trajIndices1(trajInd1),indexStrain)= cellStrain(:,2);
        end
    end
end
   
figure%('visible','off');
imshow(stateMatrix1+1,[1,1,1;COLORS_STATES])
%%%%%%%%%%%%%%%%%
createStatistics(cellClusterM1,indexStates,indexCC_time,numStates,COLORS_STATES,mu);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
table1=dataset(cellCoordinates(:,4),cellCoordinates(:,1),cellCoordinates(:,2),cellCoordinates(:,3),cellClusterM1(:,indexStates),trajIndex1,cellClusterM1(:,indexStrain));
table1.Properties.VarNames={'Image_Number','X_local','Y_local','Number_object_Number','Cluster','Traj','Strain'};
%%%%%%%%%%%%%%%%%%%%%%%%LABEL NUCLEI with their HMM STATES%%%%%%%%%%%%%%%%%    
[stateNuclei,probNuclei]=createNucleiStatePlots(m1Data,COLORS_STATES,table1,numStates,[1],pathToMovie,1);


