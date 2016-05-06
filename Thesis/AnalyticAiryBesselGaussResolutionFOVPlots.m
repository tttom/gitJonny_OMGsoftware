%%% Analytical Gaussian, Bssel, and Airy Resolution/FOV Plots

lambda=0.532; %um
n=1.33;
%resolution and FOV (NA,alpha,beta)
rG=@(NA) lambda/2./NA/0.88;
FovG=@(NA) 4*lambda*n./(NA.^2);
rG=@(NA) lambda/2./NA/0.88;
rB=@(NA,b) 0.05*pi*lambda./NA./b;
FovB=@(NA,b) lambda/2/n./b./(1-sqrt(1-(NA/n).^2));
rA=@(NA,a) max(0.12.*a*lambda/2./NA,lambda/2./NA/0.88);
FovA=@(NA,a) 6.*a*lambda/n./(1-sqrt(1-(NA/n).^2));

%resolution as function of FOV[fov],alpha,beta)
rG_FOV=@(fov) lambda/0.88/2.*sqrt(fov/4/lambda/n);
rB_FOV=@(fov,b) 0.05*pi*lambda./b/n./sqrt(1-(1-(lambda/2/n./b./fov)).^2);
rA_FOV=@(fov,a) min(lambda.*a*0.12/n./sqrt(1-(1-(6.*a*lambda/n./fov)).^2),lambda/0.88/2.*sqrt(fov/4/lambda));

% rG(0.42)
% FovG(0.42)
% rB(0.42,[0.1,0.05])
% FovB(0.42,[0.1,0.05])
% rA(0.42,[2.5,5,10,20])
% FovA(0.42,[2.5,5,10,20])

NAs=[0.01:0.01:1.33];

figure;semilogy(NAs,rG(NAs),NAs,rB(NAs,0.1),NAs,rB(NAs,0.05),NAs,rA(NAs,10),NAs,rA(NAs,20));
title('resolution vs NA');
xlabel('NA');
ylabel('Axial resolution [um]');
axis square;
xlim([0 1.33]);
ylim([10^-1 10^2]);
legend('Gaussian','Bessel10','Bessel5','Airy10','Airy20');

figure;semilogy(NAs,FovG(NAs),NAs,FovB(NAs,0.1),NAs,FovB(NAs,0.05),NAs,FovA(NAs,10),NAs,FovA(NAs,20));
title('FOV vs NA');
xlabel('NA');
ylabel('FOV [um]');
axis square;
xlim([0 1.33]);
ylim([10^0 10^3]);
legend('Gaussian','Bessel10','Bessel5','Airy10','Airy20');

fovs=[1:1000];

figure;plot(fovs,rG_FOV(fovs),fovs,rB_FOV(fovs,0.1),fovs(5:end),rB_FOV(fovs(5:end),0.05),fovs(25:end),rA_FOV(fovs(25:end),10),fovs(50:end),rA_FOV(fovs(50:end),20));
title('resolution vs FOV');
xlabel('FOV [um]');
ylabel('Axial resolution [um]');
axis square;
xlim([0 1000]);
ylim([0 10]);
legend('Gaussian','Bessel10','Bessel5','Airy10','Airy20');



