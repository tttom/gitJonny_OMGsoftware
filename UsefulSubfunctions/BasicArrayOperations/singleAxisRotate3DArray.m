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
        inputArray=zeros(9,9,9);
        inputArray(:,5,5)=1;
        inputArray(5,:,5)=1;
        inputArray(5,5,:)=1;
    end
    if nargin<2
        %default rotation axis (dimension)
        rotAxis=3;
    end
    if nargin<3
        % deault rotation angle
        rotAngle=45/360*2*pi; %radians
    end
    
    %permute





end