%%% Function:           rotate2DArray
%%% Author:             Jonathan Nylk
%%% Created:            15/07/2015
%%% Description:        This function rotates a 2-dimensional array around
%%%                     its central axis by a specified angle (in radians).
%%%
%%% Updates (latest first):
%%%
%%%
%%% END

function [outputArray]=rotate2DArray(inputArray,rotAngle)

    % initialise default values
    if nargin<1
        % default array if no array specified (after debugging change to
        % error message for no array)
        inputArray=zeros(99,99);
        inputArray(:,50)=1;
        inputArray(50,:)=1;
    end
    if nargin<2
        % deault rotation angle
        rotAngle=10/360*2*pi; %radians
    end
    
    %generate zero-centred coordinate system
    [rowRange,colRange]=meshgrid(((1:size(inputArray,2))-floor(size(inputArray,2)/2)),((1:size(inputArray,1))-floor(size(inputArray,1)/2)));
    
    %generate rotated coordinates
    rotRowRange=rowRange.*cos(rotAngle)+colRange.*sin(rotAngle);
    rotColRange=-1*rowRange.*sin(rotAngle)+colRange.*cos(rotAngle);

    %perform sub-pixel interpolated rotation
    outputArray=interp2(rowRange,colRange,inputArray,rotRowRange,rotColRange,'*linear*',0);

end