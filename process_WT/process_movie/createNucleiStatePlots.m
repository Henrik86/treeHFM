function [ stateNuclei,probNuclei ] = createNucleiStatePlots(mData,COLORS_STATES,table,numStates,indexSite,pathToMovie,resizeFactor)

stateNuclei={};
probNuclei={};
for i=1:numStates
       stateNuclei{i}= {};
       probNuclei{i}=[];
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(length(indexSite)>1)%if there are multipe sites
    numFrames=indexSite(2)-1;
else
    numFrames=length(mData.handles.Measurements.Image.Count_Nuclei);
end
numFrames=length(mData.handles.Measurements.Image.Count_Nuclei);
group_number=[];
image_number=[];
for iS=1:numFrames
      [numNuclei,c]=size(mData.handles.Measurements.Nuclei.Number_Object_Number{iS});
      group_number=[group_number;repmat(mData.handles.Measurements.Image.Group_Number{iS},numNuclei,1)];
      image_number=[image_number;repmat(mData.handles.Measurements.Image.Group_Index{iS},numNuclei,1)];
end
image_number=double(image_number);


cellstrObjectNames = fieldnames(mData.handles.Measurements.Image);
matOkObjectIX = find(~cellfun(@isempty,strfind(cellstrObjectNames,'FileName')) ,1,'first');
strObjectName = cellstrObjectNames{matOkObjectIX};
allFileNames=mData.handles.Measurements.Image.(strObjectName)';

COLORS_STATES_PLOT=[0,0,0;COLORS_STATES];
if(length(indexSite)>1)%if there are multipe sites
    numFrames=indexSite(2)-1;
else
    numFrames=length(mData.handles.Measurements.Image.Count_Nuclei)
end
numSites=length(indexSite);

se = strel('disk',5);
for iFrame=1:numFrames
     tableFrame=table(table.Image_Number==iFrame,:);
     for iSite=1:numSites
            strcat('iFrame- ',num2str(iFrame));
            traj=tableFrame.Traj;
            colorsCells=COLORS_STATES_PLOT(tableFrame.Cluster+1,:);
            hmmStates=tableFrame.Cluster;
            %img=imread(strcat(path,'/',allFileNames{iFrame}));
            %strcat(pathToFile,strcat('Nuc_',num2str(iFrame),'.tiff'))
            img=imread(strcat(pathToMovie,'/',allFileNames{iFrame}));
            if(max(double(img(:)))>1)
                %imgT=log(double(img))
                img=double(img)/max(double(img(:)));
                imgT=img;
                if(length(size(imgT))>2)
                    imgT=rgb2gray(imgT);
                end
            end

            %img=im2double(imread(strcat(pathToFile,'/',strcat('Nuc_',num2str(iFrame),'.tiff'))));
            %img=rgb2gray(img);
            [m n]=size(imgT);
            rgb=zeros(m,n,3);
            size(imgT)
            rgb(:,:,1)=imgT;
            rgb(:,:,2)=rgb(:,:,1);
            rgb(:,:,3)=rgb(:,:,1);
            rgbImage=rgb;%/255;
            imgR=rgbImage(:,:,1);
            imgG=rgbImage(:,:,1);
            imgB=rgbImage(:,:,1);
            level = graythresh(img);
            %BW = im2bw(img,mData.handles.Measurements.Image.Threshold_FinalThreshold_Nuclei{iFrame});
            BW = im2bw(img,level);
            BW=imopen(BW,se);
            L=bwlabel(BW);
            outlinedImage=bwperim(L);
            outlinedObjects=L;
            outlinedObjects(outlinedImage==0)=0;
            outlinedImage=L(outlinedImage==1);
            xValues=round(tableFrame.X_local(tableFrame.Image_Number==iFrame));
            yValues=round(tableFrame.Y_local(tableFrame.Image_Number==iFrame));
            idx = sub2ind(size(L), yValues,xValues );
            objectsImage=L(idx);
            for obj=1:nrows(traj)
                if(objectsImage(obj) ~= 0)
                    if(hmmStates(obj) ~= 0)
                        objectsImage(obj);
                        indexObj=find(outlinedObjects==objectsImage(obj));
                        color=colorsCells(obj,:);
                        imgR(indexObj)= color(1);
                        imgG(indexObj)=color(2);
                        imgB(indexObj)=color(3);


                        [r,c] = find(outlinedObjects==objectsImage(obj));
                        LSmall=L(min(r):max(r),min(c):max(c));
                        imgRSmall=imgR(min(r):max(r),min(c):max(c));
                        imgRSmall(LSmall==0)=1;
                        imgGSmall=imgG(min(r):max(r),min(c):max(c));
                        imgGSmall(LSmall==0)=1;
                        imgBSmall=imgB(min(r):max(r),min(c):max(c));
                        imgBSmall(LSmall==0)=1;
                        rgbImage=cat(3,imgRSmall,imgGSmall,imgBSmall);
                        rgbImage = imresize(rgbImage,resizeFactor);
                        stateNuclei{hmmStates(obj)}=[stateNuclei{hmmStates(obj)},rgbImage];
                    end
                end

            end
            
      
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
sortStateNucleis(stateNuclei,numStates)

