% [varargout]=calcRanges(rangeLengths,samplePitches,centerOffsets)
%
% returns uniformely spaced ranges of length rangeLengths(idx) with a elements spaced by
% samplePitches(idx) and centered on centerOffsets(idx). The center element
% is defined as the one in the center for an odd number of elements and the
% next one for an even number of elements. If a scalar is specified as sample pitch
% or center offset, it is used for all ranges. The default sample pitch is 1
% and the default center offset is 0.
%
% Example:
%    [xRange,yRange]=calcRanges([128 128],[1 1]*1e-6);
%
function [varargout]=calcRanges(rangeLengths,samplePitches,centerOffsets)
    nbDims=numel(rangeLengths);
    if nargin<2 || isempty(samplePitches),
        samplePitches=ones(1,nbDims);
    end
    if nargin<3 || isempty(centerOffsets),
        centerOffsets=zeros(1,nbDims);
    end
    % Make sure the vectors are of the same length
    samplePitches(end+1:nbDims)=samplePitches(end);
    centerOffsets(end+1:nbDims)=centerOffsets(end);
    
    for dimIdx=1:nbDims,
        varargout{dimIdx}=centerOffsets(dimIdx)+samplePitches(dimIdx)*([1:rangeLengths(dimIdx)]-1-floor(rangeLengths(dimIdx)/2));
    end
end