%%% Function:           determine2DCentreOfMass
%%% Author:             Jonathan Nylk
%%% Created:            21/07/2015
%%% Description:        This function determines x-, and y- coordinates for
%%%                     the centre of mass in a 2D array and outputs the
%%%                     coorinates.
%%%
%%% Updates (latest first):
%%%
%%%
%%% END

function [xCentre,yCentre]=determine2DCentreOfMass(inputArray,xRange,yRange,method)
    
    %default values
    if nargin<1
        % default array if no array specified (after debugging change to
        % error message for no array)
        inputArray=zeros(99,99);
        inputArray(50,50)=1;
    end
    if nargin<2
        xRange=1:size(inputArray,2);
    end
    if nargin<3
        yRange=1:size(inputArray,1);
    end
    if nargin<4
        method='com';
    end
    
    %ensure all array elements are positive
    inputArray=inputArray.*(inputArray>0);
    
    if strcmp(method,'com')
        %create 2D coordinate grid
        [X,Y]=meshgrid(xRange,yRange);
        %calculate CoM
        xCentre=sum(inputArray(:).*X(:))/sum(inputArray(:));
        yCentre=sum(inputArray(:).*Y(:))/sum(inputArray(:));
    elseif strcmp(method,'max')
        %determine maximum pixel value
        [~,maxPos]=max(inputArray(:));
        [yCentre,xCentre]=ind2sub(size(inputArray),maxPos);
    end

end