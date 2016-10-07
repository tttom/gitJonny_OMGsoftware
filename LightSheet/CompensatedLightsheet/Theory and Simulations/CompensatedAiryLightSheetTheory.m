%%% name:           CompensatedAiryLightSheetTheory
%%% author:         Jonathan Nylk
%%% date created:   19/07/2016
%%% description:    This function produces equations and plots based on the
%%%                 theoretical parameters for a compensated Airy
%%%                 light-sheet.
%%%
%%% updates (latest first):
%%%                 
%%%
%%%
%%% END %%%


% "real-space" Airy parameter [metres]
x0 = @(lambda,NA,alpha) nthroot(6 * pi .* alpha,3) .* lambda / 2 ./ NA;

% "real-space" compensation parameter [metres]
b0 = @(lambda,NA,sigma) sigma .* lambda / 2 ./ NA;

% compensation (using "real-space" variables) for checking Kai_sigma_alpha
% [m^-1]
Kai_b0_x0_conversion = @(lambda,NA,n,sigma,alpha) b0(lambda,NA,sigma) .* lambda / 2 / pi ./ n ./ (x0(lambda,NA,alpha).^3);

% compensation (using "real-space" variables) [m^-1]
Kai_b0_x0 = @(lambda,n,b0,x0) b0 .* lambda / 2 / pi / n ./ (x0.^3);

% compensation (using normalised variables) [m^-1]
Kai_sigma_alpha = @(lambda,NA,n,sigma,alpha) sigma .* (NA.^2) / 3 / (pi^2) ./ n ./ alpha ./ lambda;


% include "fudge-factor" of 10 to make theory consistent with simulation
Kai_b0_x0 = @(lambda,n,b0,x0) Kai_b0_x0(lambda,n,b0,x0) * 10; % [m^-1]
Kai_sigma_alpha = @(lambda,NA,n,sigma,alpha) Kai_sigma_alpha(lambda,NA,n,sigma,alpha) * 10; % [m^-1]


% scale to cm^-1
Kai_b0_x0 = @(lambda,n,b0,x0) Kai_b0_x0(lambda,n,b0,x0) / 100; % [cm^-1]
Kai_sigma_alpha = @(lambda,NA,n,sigma,alpha) Kai_sigma_alpha(lambda,NA,n,sigma,alpha) / 100; % [cm^-1]

% plot some parameters
NA = 0.42;
n = 1.33;
lambda = 532e-9;
[alpha,sigma] = meshgrid([1:0.01:20],[0:0.001:0.65]);

figure();
plot([0:0.001:0.65],Kai_sigma_alpha(lambda,NA,n,[0:0.001:0.65],3) .* ([0:0.001:0.65] .* (3.^0.0635) <= 0.6125)...
    ,[0:0.001:0.65],Kai_sigma_alpha(lambda,NA,n,[0:0.001:0.65],5) .* ([0:0.001:0.65] .* (5.^0.0635) <= 0.6125)...
    ,[0:0.001:0.65],Kai_sigma_alpha(lambda,NA,n,[0:0.001:0.65],7) .* ([0:0.001:0.65] .* (7.^0.0635) <= 0.6125)...
    ,[0:0.001:0.65],Kai_sigma_alpha(lambda,NA,n,[0:0.001:0.65],9) .* ([0:0.001:0.65] .* (9.^0.0635) <= 0.6125));
xlabel('sigma');ylabel('compnesation [cm^-^1]');
axis square;
title('Compensation (alpha=[3,5,7,9], sigma; n=1.33, lambda=532nm) [cm^-^1]');
legend('alpha = 3','alpha = 5','alpha = 7','alpha = 9','Location','NorthWest');

maximum_sigma_limit = (sigma .* (alpha.^0.0635) <= 0.6125);

figure();
contourf(alpha,sigma,maximum_sigma_limit,1);
axis square;
colormap gray;
xlabel('alpha');ylabel('sigma');
title('Compnesation limit (alpha, sigma; n=1.33, lambda=532nm) [cm^-^1]');

Kai_sigma_alpha_withLimits = Kai_sigma_alpha(lambda,NA,n,sigma,alpha) .* maximum_sigma_limit;
Kai_sigma_alpha_withLimits(Kai_sigma_alpha_withLimits == 0) = NaN;
figure();
[C,h] = contour(alpha,sigma,Kai_sigma_alpha_withLimits...
    ,[5,10,20,30,40,50,75,100,150,200,250]);
axis square;
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
set(h,'LineWidth',3);
colormap cool;
xlabel('alpha');ylabel('sigma');
title('Compensation with normalised parameters (alpha, sigma; NA=0.42, n=1.33, lambda=532nm) [cm^-^1] - Limits shown');

figure();
[C,h] = contour(alpha,sigma,Kai_sigma_alpha(lambda,NA,n,sigma,alpha)...
    ,[5,10,20,30,40,50,75,100,150,200,250]);
axis square;
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
set(h,'LineWidth',3);
colormap cool;
xlabel('alpha');ylabel('sigma');
title('Compensation with normalised parameters (alpha, sigma; NA=0.42, n=1.33, lambda=532nm) [cm^-^1]');

[x0,b0] = meshgrid([1.5:0.001:4.5]*1e-6,[0:0.002:4]*1e-7);

maximum_b0_limit = (b0 .* (x0.^(3 * 0.0635)) <= 0.6125 * (6 * pi).^0.0635 * (lambda / 2 / NA)^(3 * 0.0635 + 1));

Kai_b0_x0_withLimits = Kai_b0_x0(lambda,n,b0,x0) .* maximum_b0_limit;
Kai_b0_x0_withLimits(Kai_b0_x0_withLimits == 0) = NaN;

figure();
contourf(x0*1e6,b0*1e6,maximum_b0_limit,1);
axis square;
colormap gray;
xlabel('x0 [um]');ylabel('b0 [um]');
title('Compnesation limit (x0, b0; n=1.33, lambda=532nm) [cm^-^1]');

figure();
[C,h] = contour(x0*1e6,b0*1e6,Kai_b0_x0_withLimits...
    ,[5,10,20,30,40,50,75,100,150,200,250]);
axis square;
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
set(h,'LineWidth',3);
colormap cool;
xlabel('x0 [um]');ylabel('b0 [um]');
title('Compensation with real-space parameters (x0, b0; n=1.33, lambda=532nm) [cm^-^1] - with limits');

figure();
[C,h] = contour(x0*1e6,b0*1e6,Kai_b0_x0(lambda,n,b0,x0)...
    ,[5,10,20,30,40,50,75,100,150,200,250]);
axis square;
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
set(h,'LineWidth',3);
colormap cool;
xlabel('x0 [um]');ylabel('b0 [um]');
title('Compensation with real-space parameters (x0, b0; n=1.33, lambda=532nm) [cm^-^1]');

% plot some parameters
NA = 0.42;
n = 1.33;
lambda = 488e-9;
[alpha,sigma] = meshgrid([1:0.01:20],[0:0.001:0.65]);

maximum_sigma_limit = (sigma .* (alpha.^0.0635) <= 0.6125);

figure();
contourf(alpha,sigma,maximum_sigma_limit,1);
axis square;
colormap gray;
xlabel('alpha');ylabel('sigma');
title('Compnesation limit (alpha, sigma; n=1.33, lambda=488nm) [cm^-^1]');

Kai_sigma_alpha_withLimits = Kai_sigma_alpha(lambda,NA,n,sigma,alpha) .* maximum_sigma_limit;
Kai_sigma_alpha_withLimits(Kai_sigma_alpha_withLimits == 0) = NaN;
figure();
[C,h] = contour(alpha,sigma,Kai_sigma_alpha_withLimits...
    ,[5,10,20,30,40,50,75,100,150,200,250]);
axis square;
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
set(h,'LineWidth',3);
colormap cool;
xlabel('alpha');ylabel('sigma');
title('Compensation with normalised parameters (alpha, sigma; NA=0.42, n=1.33, lambda=488nm) [cm^-^1] - Limits shown');

figure();
[C,h] = contour(alpha,sigma,Kai_sigma_alpha(lambda,NA,n,sigma,alpha)...
    ,[5,10,20,30,40,50,75,100,150,200,250]);
axis square;
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
set(h,'LineWidth',3);
colormap cool;
xlabel('alpha');ylabel('sigma');
title('Compensation with normalised parameters (alpha, sigma; NA=0.42, n=1.33, lambda=488nm) [cm^-^1]');

[x0,b0] = meshgrid([1.5:0.001:4.5]*1e-6,[0:0.002:4]*1e-7);

maximum_b0_limit = (b0 .* (x0.^(3 * 0.0635)) <= 0.6125 * (6 * pi).^0.0635 * (lambda / 2 / NA)^(3 * 0.0635 + 1));

Kai_b0_x0_withLimits = Kai_b0_x0(lambda,n,b0,x0) .* maximum_b0_limit;
Kai_b0_x0_withLimits(Kai_b0_x0_withLimits == 0) = NaN;

figure();
contourf(x0*1e6,b0*1e6,maximum_b0_limit,1);
axis square;
colormap gray;
xlabel('x0 [um]');ylabel('b0 [um]');
title('Compnesation limit (x0, b0; n=1.33, lambda=488nm) [cm^-^1]');

figure();
[C,h] = contour(x0*1e6,b0*1e6,Kai_b0_x0_withLimits...
    ,[5,10,20,30,40,50,75,100,150,200,250]);
axis square;
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
set(h,'LineWidth',3);
colormap cool;
xlabel('x0 [um]');ylabel('b0 [um]');
title('Compensation with real-space parameters (x0, b0; n=1.33, lambda=488nm) [cm^-^1] - with limits');

figure();
[C,h] = contour(x0*1e6,b0*1e6,Kai_b0_x0(lambda,n,b0,x0)...
    ,[5,10,20,30,40,50,75,100,150,200,250]);
axis square;
set(h,'ShowText','on','TextStep',get(h,'LevelStep'));
set(h,'LineWidth',3);
colormap cool;
xlabel('x0 [um]');ylabel('b0 [um]');
title('Compensation with real-space parameters (x0, b0; n=1.33, lambda=488nm) [cm^-^1]');

% colour map
figure();imagesc(repmat([5:1:250].',[1,10]));colormap cool;

