function [ output_args ] = sortStateNucleis(allStateNuclei,numStates,outputFolderState)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
figure
maxCellNumber=50;
valueToBreak=10;
numCellsToProcesss=10;
for s=1:numStates
    subplot(numStates/2,numStates/2,s)
    %
    maxRow=0;
    maxCol=0;
    stateImages=allStateNuclei{s};
    if(~isempty(stateImages))
        stateImages;
        %find maximum number of rows, cols
        for img=1:length(stateImages)
            if(nrows(stateImages{img})>maxRow)
                maxRow=nrows(stateImages{img});
            end
            if(ncols(stateImages{img})>maxCol)
                maxCol=ncols(stateImages{img});
            end
        end
    
        numImg=length(allStateNuclei{s});
        numRows=round(numCellsToProcesss/valueToBreak);
        stateMatrix=cat(3,ones(valueToBreak*maxRow,(maxCellNumber/10)*maxCol),ones(valueToBreak*maxRow,(maxCellNumber/10)*maxCol),ones(valueToBreak*maxRow,(maxCellNumber/10)*maxCol));
        rowCounter=1;
        existingMatrix=cat(3,[1],[1],[1]);
    
        hmmStateRows={};
        for nr=1:numRows
            hmmStateRows{nr}=existingMatrix;
        end
        cellCounter=0;
        numCellsToProcesssTemp=numCellsToProcesss;
        if(numCellsToProcesssTemp>numImg)
            numCellsToProcesssTemp=numImg;
        end
        imagesToChoose=randi(numel(stateImages),1,numCellsToProcesssTemp);
        for r=1:numCellsToProcesssTemp
             strcat('State',num2str(r))
             rgbImage=stateImages{imagesToChoose(r)};

             while(nrows(existingMatrix)<nrows(rgbImage))
                 existingMatrix=[existingMatrix;cat(3,ones(1,ncols(existingMatrix)),ones(1,ncols(existingMatrix)),ones(1,ncols(existingMatrix)))];
             end
             while(nrows(existingMatrix)>nrows(rgbImage))
                 rgbImage=[rgbImage;cat(3,ones(1,ncols(rgbImage)),ones(1,ncols(rgbImage)),ones(1,ncols(rgbImage)))];
             end

            existingMatrix=[existingMatrix,cat(3,ones(nrows(rgbImage),1),ones(nrows(rgbImage),1),ones(nrows(rgbImage),1)),rgbImage];
            
            if(cellCounter>numRows)
                hmmStateRows{rowCounter}=existingMatrix;
                existingMatrix=cat(3,[1],[1],[1]);
                rowCounter=rowCounter+1;
                cellCounter=0;
                
            end
            cellCounter=cellCounter+1;
            if(r>maxCellNumber)
             %   break
            end
        end
        'write state matrix';
        stateCells=cat(3,[1],[1],[1]);
        rowMatrix=1;
        maxColMatrix=1;
        for nr=1:length(hmmStateRows)
            row=hmmStateRows{nr};
            if(~isempty(row))
                if(ncols(row)>maxColMatrix)
                    maxColMatrix=ncols(row);
                end
                stateMatrix(rowMatrix:(rowMatrix+nrows(row)-1),1:ncols(row),:)=hmmStateRows{nr};
                rowMatrix=rowMatrix+nrows(row)-1;
            end
        end
        stateMatrix=stateMatrix(:,1:maxColMatrix,:);
        stateMatrix=stateMatrix(1:rowMatrix,:,:);
    end
        subimage(stateMatrix)
        set(gca,'ytick',[]) 
        set(gca,'xtick',[])
end

end

