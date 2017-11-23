function [ output_args ] = createStatistics( cellClusterM1,indexStates,indexCC_time,numStates,COLORS_STATES,muN,stateMatrix1)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%%%%%%%
indexCell=1:nrows(cellClusterM1);
cellClusterMT=cellClusterM1;
indToDelete=find(cellClusterM1(:,indexStates)==0);
indexCellDT=indexCell;
indexCellDT(indToDelete)=[];
cellClusterMT(indToDelete,:)=[];
statesMT=cellClusterMT(:,indexStates);
colorsPlot=[]
for nC=1:nrows(cellClusterMT) 
    colorsPlot=[colorsPlot;COLORS_STATES(cellClusterMT(nC,indexStates),:)];
end
cellClusterMT=[muN';cellClusterMT(:,1:(indexStates-2))];
[pcF,scoreF,latentF,tsquareF]=princomp(cellClusterMT);
%%%%%
figure
subplot(2,2,1);
hold on
xlim([min(scoreF(:,1)) max(scoreF(:,1))])
scatter(scoreF(7:nrows(cellClusterMT),1),scoreF(7:nrows(cellClusterMT),2),3,colorsPlot,'+');
scores1=scoreF(7:nrows(cellClusterMT),1);
scores2=scoreF(7:nrows(cellClusterMT),2);
for nH=1:numStates
    [ x_bin_centers, y_bin_centers,counts ] = createHistMatrix( scores1(statesMT==nH), scores2(statesMT==nH),10 );
    [C H] = contour(x_bin_centers, y_bin_centers, counts',2,'Linecolor',COLORS_STATES(nH,:));
    set (H, 'LineWidth', 2); 
end
scatter(scoreF(1:6,1),scoreF(1:1:6,2),100,COLORS_STATES(1:6,:),'o','filled','MarkerEdgeColor','k');
%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT HISTOGRAMS%%%%%%%%%%%%%%%%%%%%%%%%%%
fontSize=12;
%%%
subplot(2,2,2);
nHistStates=10;
indicesSwitched=[(nHistStates/2):nHistStates,1:((nHistStates/2)-1)];
    hold on
    for indS1=1:numStates
        indS1;
        timepointsState=cellClusterM1(cellClusterM1(:,indexStates)==indS1,indexCC_time);
        [h,b]=hist(timepointsState,[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]);
        h=h(indicesSwitched);
        b=b(indicesSwitched);
        p=plot(h,'-o','LineWidth',1,...
                'MarkerEdgeColor',COLORS_STATES(indS1,:),...
                'MarkerFaceColor',COLORS_STATES(indS1,:),...
                'MarkerSize',5);
         set(gca,'XTickLabel', b);
       set(p,'Color',COLORS_STATES(indS1,:),'LineWidth',4)
       set(gca,'FontSize',fontSize)

    end
    hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%STATE NUMBERS%%%%%%%%%%%%%%%%%%%%%%%%%%
stateNumber=tabulate(cellClusterM1(:,indexStates));
stateNumber(stateNumber(:,1)==0,:)=[];
allStates1=sum(stateNumber(:,2));
subplot(2,2,3);
hold on
maxBarValue=max(stateNumber(:,2)/allStates1);
counter=1;
for n=1:nrows(stateNumber)
    b1=bar(counter,stateNumber(n,2)/allStates1,'FaceColor',COLORS_STATES(stateNumber(n,1),:))
    counter=counter+1;
end
hold off
set(gca,'xTickLabel',[' '])
set(gca,'FontSize',fontSize)
end