function [zAmp,zOffset,xCentre,yCentre,xWidth,yWidth] = GaussianFit2D(input_matrix,debug_flag)
    % Fits a Gaussian function to the data within the array
    % Relies on getting a 2d array as input matrix
    
    if nargin > 2
        debug_flag = 0;
    end
    
    % Define coordiante system (in pixels)
    dim1Range = [1:size(input_matrix,1)];
    dim2Range = [1:size(input_matrix,2)];
    [dim1Coords,dim2Coords] = ndgrid(dim1Range,dim2Range);
    
    % Define "D Gaussian function
    % params = [(1) amplitude, (2) dim1_centre, (3) dim2_centre, (4) dim1_width, (5) dim2_width, (6) y_offset];
    functor_2DGaussian = @(dim1Coords,dim2Coords,params) params(1) ...
        .* exp(-1 .* ((dim1Coords - params(2)).^2) ./ 2 ./ (params(4)).^2) ...
        .* exp(-1 .* ((dim2Coords - params(3)).^2) ./ 2 ./ (params(5)).^2) ...
        + params(6);
    
    % Function to minimise to fit 2D Gaussian to data
    functor_2DGaussian_minimisation = @(dim1Coords,dim2Coords,params,data)...
        sum(sum((functor_2DGaussian(dim1Coords,dim2Coords,params) - data).^2));
    
    % input guess params
    params0 = zeros([1,6]);
    [params0(1),max_Idx] = max(input_matrix(:)); % guess Gaussian amplitude
    [params0(2),params0(3)] = ind2sub(size(input_matrix),max_Idx); % guess of dim1/2 centres
    params0(4:5) = 2; % guess of dim1/2 widths
    params0(6) = min(input_matrix(:)); % guess z_offset
    
    % fitting
    [params,~,exitflag] = fminsearch(@(params) functor_2DGaussian_minimisation(dim1Coords,dim2Coords,params,input_matrix),params0);
    if exitflag ==1
        % if converged on a solution, use these values
    else
        % if can't fit, use initial guesses
        params = params0;
    end
    zAmp = params(1);
    yCentre = params(2);
    xCentre = params(3);
    yWidth = params(4);
    xWidth = params(5);
    zOffset = params(6);
    
    if debug_flag
        figure;
        subplot(2,2,1);
        imagesc(dim1Range,dim2Range,input_matrix);
        axis image;
        title('input matrix');
        subplot(2,2,2);
        imagesc(dim1Range,dim2Range,functor_2DGaussian(dim1Coords,dim2Coords,params));
        axis image;
        title('fitted 2D Gaussian');
        subplot(2,2,3);
        plot(dim2Range.',input_matrix(round(params(2)),:).',dim2Range.',functor_2DGaussian(round(params(2)),dim2Coords,params).');
        title('x cross-section');
        subplot(2,2,4);
        plot(dim1Range,input_matrix(:,round(params(3))),dim1Range,functor_2DGaussian(dim1Coords,round(params(3)),params));
        title('y cross-section');
    
end


