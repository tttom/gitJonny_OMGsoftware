%%% Function:           rotate2DArray
%%% Author:             Jonathan Nylk
%%% Created:            15/07/2015
%%% Description:        This function rotates a 2-dimensional array around
%%%                     its central axis by a specified angle (in radians)
%%%                     anti-clockwise.
%%%
%%% Updates (latest first):
%%% 19/20/2015:         If a zero-centred coordinate system is given as an
%%%                     input, this is rotated and output.
%%%
%%% END

function [outputArray,outputColRange,outputRowRange]=rotate2DArray(inputArray,rotAngle,inputColRange,inputRowRange,cubic_FLAG,displayRotation)

    % initialise default values
    if nargin<1
        % default array if no array specified (after debugging change to
        % error message for no array)
        inputArray=zeros(200,100);
        inputArray(:,50)=1;
        inputArray(100,:)=1;
    end
    if nargin<2
        % deault rotation angle
        rotAngle=45/360*2*pi; %radians
    end
    
    if nargin<4
        %generate zero-centred coordinate system
        [rowRange,colRange]=meshgrid(((1:size(inputArray,2))-floor(size(inputArray,2)/2)),((1:size(inputArray,1))-floor(size(inputArray,1)/2)));
    else
        [rowRange,colRange]=meshgrid(inputRowRange,inputColRange);
        clear inputRowRange inputColRange;
    end
    
    if nargin<5
        % default linear interpolation
       cubic_FLAG=0; 
    end
    
    if nargin<6
        % default display rotated image
       displayRotation=1; 
    end
    
    %check rowRange and colRange are zero centred, correct if not
    if rowRange(1,floor(size(rowRange,2)/2))~=0
        rowOffset=rowRange(1,floor(size(rowRange,2)/2));
        rowRange=rowRange-rowOffset;
    else
        rowOffset=0;
    end
    if colRange(floor(size(colRange,1)/2),1)~=0
        colOffset=colRange(floor(size(colRange,1)/2),1);
        colRange=colRange-colOffset;
    else
        colOffset=0;
    end
    
    % row and column interpixel spacing
    rowStepSize=rowRange(1,2)-rowRange(1,1);
    colStepSize=colRange(2,1)-colRange(1,1);
    
    %generate rotated coordinates
    rotRowRange=rowRange.*cos(rotAngle)+colRange.*sin(rotAngle);
    rotColRange=-1*rowRange.*sin(rotAngle)+colRange.*cos(rotAngle);

    %rotate coordinate offset
    rotRowOffset=rowOffset.*cos(rotAngle)+colOffset.*sin(rotAngle);
    rotColOffset=-1*rowOffset.*sin(rotAngle)+colOffset.*cos(rotAngle);

    %add offset back into coordinates
    rotRowRange=rotRowRange+rotRowOffset;
    rotColRange=rotColRange+rotColOffset;
    
    %perform sub-pixel interpolated rotation
    if cubic_FLAG==1
        outputArray=interp2(rowRange,colRange,inputArray,rotRowRange,rotColRange,'*cubic*',0);
    else
        outputArray=interp2(rowRange,colRange,inputArray,rotRowRange,rotColRange,'*linear*',0);
    end
    
    %output coords
    outputColRange=rotColRange(:,floor(size(rotColRange,2)/2));
    outputRowRange=rotRowRange(floor(size(rotRowRange,1)/2),:);
    
    if displayRotation
        figure();
        subplot(1,2,1);imagesc(rowRange(floor(size(rowRange,1)/2),:),colRange(:,floor(size(colRange,2)/2)),inputArray);axis image;
        subplot(1,2,2);imagesc(rotRowRange(floor(size(rotRowRange,1)/2),:),rotColRange(:,floor(size(rotColRange,2)/2)),outputArray);axis image;
    end
    
end