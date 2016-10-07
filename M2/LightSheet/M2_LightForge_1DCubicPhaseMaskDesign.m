function M2_LightForge_1DCubicPhaseMaskDesign()

% Produces the design of a 1D cubic phase mask with specified alpha value
% over specified x,y range, and outputs the data in the required GridXYZ
% format.

    fullWidth=7.5e3; %in um
    halfWidth=fullWidth/2;

    x=[-halfWidth:10:halfWidth];
    y=[-halfWidth:10:halfWidth];
    
    normalised_x=x/halfWidth;
    normalised_y=y/halfWidth;
    
    [normX,normY]=meshgrid(normalised_x,normalised_y);
    
%     alpha=6.1447;
    alpha=10.7533;
    cubic_profile=alpha.*(normX.^3);
    
    %make all positive values
    cubic_profile=cubic_profile-min(cubic_profile(:));
    
    %prepare in matrix form:
    gridxyz=zeros(size(cubic_profile)+1,'double');
    
    gridxyz(1,2:end)=x;
    gridxyz(2:end,1)=y;
    gridxyz(2:end,2:end)=cubic_profile;
    
%     imagesc(gridxyz);axis image;
    imagesc(cubic_profile);axis image;
    
    dlmwrite('gridxyz.txt',gridxyz,'delimiter','\t','precision','%.3f','newline','pc');
    
end
