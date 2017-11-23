function [ output_args ] = plotTransitionMatrices( transProbSeqN,transProbDivN,COLORS_STATES,indicesXD )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%%%%%%%%%COLORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hotColors=hot(100);
colors_hot=hotColors(100:-1:30,:);
fontSize=12;
numStates=nrows(transProbSeqN);
%%%%%%%%%%% write transition matrix
%%%%%%%%%%% Sequences%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
figure
transmat1_nd=transProbSeqN - diag(diag(transProbSeqN));

g2=subplot(2,2,(1:2));
subimage(1:numStates,COLORS_STATES)
colormap(COLORS_STATES)
set(gca,'xtick',[],'ytick',[])
p = get(g2,'position');
p(4) = p(4)*0.35;
p(3) = 0.2417;% Add 10 percent to height
p(2) = p(2)-0.13;
p(1) = 0.5703;
set(g2, 'position', p);
set(gca,'xtick',[],'ytick',[])

g3=subplot(2,2,3);
subimage((1:numStates)',COLORS_STATES)

p = get(g3,'position');
p(1) = p(1)+0.27;
p(3) = p(3)*0.5; % Add 10 percent to height
set(g3, 'position', p);
set(gca,'xtick',[],'ytick',[])
g4=subplot(2,2,4);
colormap(colors_hot)
set(gca,'xtick',[],'ytick',[])
diffMatrix=abs(transmat1_nd);
imagesc(diffMatrix)
colorbar('location','EastOutside')
p = get(g4,'position');
set(gca,'xtick',[],'ytick',[])
set(gca,'FontSize',fontSize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
transmat1_nd=transProbSeqN;
D = diag(diag(transmat1_nd));
diagVector=diag(transmat1_nd);
diagM=repmat(diagVector,1,ncols(transmat1_nd));
eBii=1./(1-D);
pBij=transmat1_nd./(1-diagM);
pBij=pBij-diag(diag(pBij));
matrixAll=eBii;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% write transition matrix Division%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,1,2);
transmat1_nd=transProbDivN;
%%%%%
numStatesDiv=numStates*numStates;
%%%%%
g2=subplot(2,2,(1:2));
colors=COLORS_STATES(1:numStates,:);
colorsLeft=repmat(colors,numStates,1);
colorsLeft=[];
for i = 1 : numStates
    colorsLeft=[colorsLeft;repmat(COLORS_STATES(i,:),numStates,1)];
end
colorsRight=colorsLeft(indicesXD,:);
colorsAll=[colorsLeft;colorsRight];
subimage([1:numStatesDiv;(numStatesDiv+1):(numStatesDiv*2)],colorsAll)

colormap(repmat(colorsAll,3,1))
set(gca,'xtick',[],'ytick',[])
p = get(g2,'position');
p(4) = p(4)*0.35;
p(3) = 0.2417;% Add 10 percent to height
p(2) = p(2)-0.13;
p(1) = 0.5703;
set(g2, 'position', p);
set(gca,'xtick',[],'ytick',[])

g3=subplot(2,2,3);
subimage((1:numStates)',COLORS_STATES)
colormap(COLORS_STATES)
p = get(g3,'position');
p(1) = p(1)+0.27;
p(3) = p(3)*0.5; % Add 10 percent to height
set(g3, 'position', p);
set(gca,'xtick',[],'ytick',[])
g4=subplot(2,2,4);
colormap(colors_hot)
set(gca,'xtick',[],'ytick',[])
diffMatrix=transmat1_nd;
imagesc(diffMatrix)
colorbar('location','EastOutside')
p = get(g4,'position');
set(gca,'xtick',[],'ytick',[])
set(gca,'FontSize',fontSize)
set(gcf, 'renderer','opengl'); % or painters, or zbuffer 
%%%%%%%%%%%%%%%%%%%%% NUM TRANSITIONS DIV %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

