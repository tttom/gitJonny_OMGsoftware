%%% Function:           singleAxisRotate3DArray
%%% Author:             Jonathan Nylk
%%% Created:            15/07/2015
%%% Description:        This function rotates a 3-dimensional array around
%%%                     the centre of one axis (dimension) by a specified 
%%%                     angle (in radians).
%%%
%%% Updates (latest first):
%%%
%%%
%%% END

function [outputArray]=singleAxisRotate3DArray(inputArray,rotAxis,rotAngle)

    % initialise default values
    if nargin<1
        % default array if no array specified (after debugging change to
        % error message for no array)
        inputArray=zeros(99,99,99);
        inputArray(:,50,50)=1;
        inputArray(50,:,50)=1;
        inputArray(50,50,:)=1;
    end
    if nargin<2
        %default rotation axis (dimension)
        rotAxis=2;
    end
    if nargin<3
        % deault rotation angle
        rotAngle=10/360*2*pi; %radians
    end
    
    %permute rotation axis to 3rd position
    outputArray=permute(inputArray,circshift([1;2;3],3-rotAxis)');
    for zIdx=1:size(outputArray,3)
        outputArray(:,:,zIdx)=rotate2DArray(outputArray(:,:,zIdx),rotAngle);
    end
    
    %inverse permute array back to original orientation
    outputArray=ipermute(outputArray,circshift([1;2;3],3-rotAxis)');

end