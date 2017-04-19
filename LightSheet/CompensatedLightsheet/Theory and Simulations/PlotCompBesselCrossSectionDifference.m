fval_x = zeros(1,151+100);
RMS_x = zeros(1,151+100);
A_x = zeros(1,151+100);
for xIdx = 151-100:151+100
    
    [A_x(xIdx),fval_x(xIdx),residual] = BesselCrossSectionDataOptimisation(BesselProfile_Sigma0_00_Sheet(:,xIdx),BesselProfile_Sigma0_11_Sheet(:,xIdx));
    
    RMS_x(xIdx) = sqrt(mean((BesselProfile_Sigma0_00_Sheet(:,xIdx) - A_x(xIdx) .* BesselProfile_Sigma0_11_Sheet(:,xIdx)).^2));
    
    
    if max(xIdx == [151-100,151-50,151,151+50,151+100])
        figure;
        subplot(1,2,1);
        plot(zRange,BesselProfile_Sigma0_00_Sheet(:,xIdx),zRange,A_x(xIdx) .* BesselProfile_Sigma0_11_Sheet(:,xIdx));
        xlabel('z [um]'); ylabel('Intensity [a.u.]');
        xlim([-20 20]);
        title(strcat('x = ',num2str(xRange(xIdx))));
        axis square;

        subplot(1,2,2);
        plot(zRange,residual);
        xlabel('z [um]'); ylabel('Intensity [a.u.]');
        xlim([-20 20]);
        axis square;
        drawnow;shg;
    end
    
end

figure;
[Ax,H1,H2] = plotyy(xRange(151-100:151+100),RMS_x(151-100:151+100),xRange(151-100:151+100),A_x(151-100:151+100));
axes(Ax(1));
xlabel('x [um]');
ylabel('RMS error [a.u.]');
axes(Ax(2));
xlabel('x [um]');
ylabel('A [a.u.]');