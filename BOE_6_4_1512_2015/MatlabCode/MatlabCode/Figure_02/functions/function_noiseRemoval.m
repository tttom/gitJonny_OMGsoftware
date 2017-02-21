function [movOut,E,PP,val0] = function_noiseRemoval(movIn, N, EigVecUsed)
% function [movOut,E,PP] = function_noiseRemoval_v2(movIn, N)
%
% PURPOSE:
%    Remove background from movie.
%
% INPUT:
%    movIn: Three-dimensional array [Image rows, Image columns, Frame
%    number].
% 
%    fv: Fraction of variability that should be included in movOut. Range 
%    is [0, 1], the default is [0.9].
% 
%    plotting: Logical true or false input determines whether the effect of 
%    this function is plotted.
% 
% OUTPUT:
%    movOut:  Three-dimensional array [Image rows, Image columns, Frame
%    number].
% 
% CREATED: Martin Verner Gammelgaard Kristensen, St Andrews, January 2014.
%
%%
if nargin < 2,	N = 10;	end
if nargin < 3,	EigVecUsed = 1:N;	end

% Parameters
height          = size(movIn,1);      %Number of rows in each frame
width           = size(movIn,2);      %Number of columns in each frame
nFrames         = size(movIn,3);      %Number of frames

vec0            = zeros(numel(movIn(:,:,1)),N);
val0            = zeros(N,N);
meanmat         = zeros(numel(movIn(:,:,1)),1);
for J = 1:floor(nFrames/50)
    frIND           = 1+(J-1)*50:J*50;
    mat             = reshape(movIn(:,:,frIND),numel(movIn(:,:,1)),[]);
    meanmat         = (meanmat*(J-1) + mean(mat,2))/J;
    mat0            = vec0*sqrt(val0);
    mat1            = mat - repmat(meanmat,1,50);
    mat2            = [mat0,mat1];
    m               = mat2'*mat2;
    [vec val]       = eigs(m,N);
    vec0 = mat2*vec; val0 = val;
    for l = 1:N
        vec0(:,l)            = vec0(:,l)/norm(vec0(:,l));
    end

    % Progress
    if J/round(floor(nFrames/50)/100) == round(J/round(floor(nFrames/50)/100))
        clc
        display([num2str(100*J/floor(nFrames/50)),'%'])
    end
end

%% Last <50 frames
if 50*floor(nFrames/50) < nFrames
    frIND           = floor(nFrames/50)+1:nFrames;
    mat             = reshape(movIn(:,:,frIND),numel(movIn(:,:,1)),[]);
    meanmat         = (meanmat*floor(nFrames/50) + mean(mat,2))/nFrames;
    mat0            = vec0*sqrt(val0);
    mat1            = mat - repmat(meanmat,1,length(frIND));
    mat2            = [mat0,mat1];
    m               = mat2'*mat2;
    [vec val]       = eigs(m,N);
    vec0 = mat2*vec; val0 = val;
    for l = 1:N
        vec0(:,l)            = vec0(:,l)/norm(vec0(:,l));
    end
end

%% Reconstructing movie
    A               = reshape(movIn,numel(movIn(:,:,1)),[]);
    ma              = mean(A');
    mA              = repmat(ma,size(A,2),1)';
    dA              = A-mA;

    B               = vec0(:,EigVecUsed)*(dA'*vec0(:,EigVecUsed))' + mA;
    movOut          = reshape(B,height,width,nFrames);
    E               = reshape(vec0,height,width,N);
    PP              = cumsum(diag(val0)/trace(m));

end