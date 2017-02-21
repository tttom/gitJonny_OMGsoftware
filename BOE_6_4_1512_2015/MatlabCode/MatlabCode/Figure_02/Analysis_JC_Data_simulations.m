% CREATED: Martin Verner Gammelgaard Kristensen, University of St Andrews, May 2014.
 
%% * Clearing memory
close all
clear
clc

%% * Making a list of individual file names in the "DataFolder".
DataFolder                  = 'Jonny and Claire - Data\2014_05_22\800nmGRINlens\3um';             %Folder containing the data files
FOLDER                      = [cd, filesep, DataFolder];                            %Full path to folder containing the data files
TYPE                        = 'mat';                                                %File format
% SIZE                        = [1, 100];                                             %Range of files size in MB
% [FN,~]                      = filenames(FOLDER,TYPE,SIZE*10^3);                     %Function which produce a list (string) of the file names

%% * Loading DATA
load([FOLDER,filesep,'DATA_F1_3um'])

%% * Parameters
pTF                         = logical([0 1 1]);                         %Parameters to fit
plotFLAG                    = true;
DATA.settings.sectionFLAG 	= 0;
DATA.settings.nSections    	= 20;

% Derived sampling parameters 
DATA.settings.df        	= 1/DATA.settings.nFrames;                  %Frequency resolution
DATA.settings.dt        	= 1/DATA.settings.fps;                      %Temporal resolution

% Setting up additional variables
scrnsz                      = [100,100,1700,1000];%get(0,'ScreenSize');
d                           = 3e-6;                                     %Bead diameter [1.90 3.00 4.17]*1e-6
sdtd                        = 2e-8;                                     %Bead diameter standard deviation [0.02 0.02 0.03]*1e-6 
r                           = d/2;                                      %Bead radius
sdtr                        = sdtd/2;                                	%Bead radius standard deviation
s                           = 250*1E-9;%[142, 190, 237, 142, 142; 427, 735, 427, 356, 213; 640, 569, 403, 735, 332]*1E-9;
stds                        = 200*1e-9;
h                           = s+r;
a1 = 9/16; a2 = 1/8; a3 = 45/256; a4 = 1/16;
A                           = (1 - a1*r./h + a2*(r./h).^3 - a3*(r./h).^4 - a4*(r./h).^5) / r;
nFrames                     = DATA.settings.nFrames;
fps                       	= DATA.settings.fps;
t                           = (0:nFrames-1)*DATA.settings.dt;           %Time vector
p2m                         = 0.1766*1e-6;                              %Pixel to meter
kB                          = 1.3806488e-23;                            %Boltzmann constant
T                           = 25.6+272.15;%DATA.settings.T;             	%Temperature
eta                         = viscosityWater(T);                        %Viscosity of water at T [Pa s]
a0                          = 6*pi*eta;
gamma                       = a0/A;
Dcalc                       = kB*T/gamma;

% Calculating standard deviations
dgammadr                    = -(a0/A^2) * (-1/r^2 + a1/h^2 + a2*r*(2*s-r)/h^4 - a3*r^2*(3*s-r)/h^5 - a4*r^3*(4*s-r)/h^6);
dgammads                    = -(a0/A^2) * (a1/h^2 - 3*a2*r^2/h^4 + 4*a3*r^3/h^5 + 5*a4*r^4/h^6);
sdtgamma                    = sqrt(dgammadr^2*sdtr^2 + dgammads^2*stds^2);

% Fitting parameters - start values
fc                          = 0;                                        %Corner frequency
E                           = 0;                                        %Level of white background noise in the power spectrum
D                           = Dcalc;%1e-13;                                    %Einstein's diffusion constant
gamma                       = kB*T/D;                                   %
fPar                        = [1 E D];                                  %Fitting parameters, [fit.fc fit.E fit.D]
W                           = 400*1e-6;                                 %Camera exposure time, [s]. If "DATA.settings.W = 0" then parameter is not used
nA                          = 3;                                        %n=-nA:nA in anti-aliasing sum of the fitmodel. Thue typically used nA = 3;
cPar                        = [fps W nA];                               %Constant parameters, [DATA.settings.F DATA.settings.fit.W DATA.settings.fit.nA];

% Fitting parameters - saving
DATA.settings.d             = d;
fitting.fc              	= fc;
fitting.E                 	= E;
fitting.D                 	= D;
fitting.fPar             	= fPar;
fitting.pTF             	= pTF;
fitting.W                  	= W;
fitting.nA              	= nA;
fitting.cPar             	= cPar;
DATA.settings.fitting       = fitting;

% Temperature and Faxen's correction estimation
TfitRough                   = D*gamma/kB-272.15;
FaxenObserved               = Dcalc/D;

%% Temporary test
k               = 0;
ini          	= 1;
for ID = ini:5:ini+25
    k = k+1;
    track{k} = [DATA.DFT(ID).r'-0.5, DATA.DFT(ID).c'+0.5];
    
end
range           = 1:10000;
range           = range + 1e3;

scrnsz          = get(0,'ScreenSize');
fig1            = figure;
set(fig1,'units','pixels','Position', scrnsz)
subplot(1,6,1),plot(track{1}(range,:)),title(num2str(std(track{1}(range,:))))
subplot(1,6,2),plot(track{2}(range,:)),title(num2str(std(track{2}(range,:))))
subplot(1,6,3),plot(track{3}(range,:)),title(num2str(std(track{3}(range,:))))
subplot(1,6,4),plot(track{4}(range,:)),title(num2str(std(track{4}(range,:))))
subplot(1,6,5),plot(track{5}(range,:)),title(num2str(std(track{5}(range,:))))
subplot(1,6,6),plot(track{6}(range,:)),title(num2str(std(track{6}(range,:))))
return
%% * For loop:
for ID = 1:size(DATA.DFT,2)
    %% -        1. Sectioning data
    %              Data is sectioned into (n) parts of equivalent length
    if DATA.settings.sectionFLAG == 1 || DATA.settings.sectionFLAG == 2
        n                               = floor(DATA.settings.nFrames/DATA.settings.nSections);
        
        for l = 1:DATA.settings.nSections
            if DATA.settings.sectionFLAG == 1 %downsample measuring time
                index                           = n*(l-1)+1:n*l;
            elseif DATA.settings.sectionFLAG == 2 %downsample sampling frequency
                index                           = l:DATA.settings.nSections:n*DATA.settings.nSections+l-10;
            end
            DATA.DFT(ID).section.r(:,l)   	= DATA.DFT(ID).r(index);
            DATA.DFT(ID).section.c(:,l)    	= DATA.DFT(ID).c(index);
            if l == 1; ts = t(index); end %Sectioned time vector
        end
    end

	%% -        2. Fourier transforming data
    %              Both full (S) and singlesided (Sss) Fourier transform, 
    %              as well as the power spectral density (PSD), is found 
    %              for both the full signals (Data{x}.s.full) and the
    %              sectioned signals (Data{x}.s.sectioned).
    prefiltering                                                                                                    = 2;
    [~, DATA.DFT(ID).R, fss, DATA.DFT(ID).Rss, DATA.DFT(ID).Rpsd]                                                   = transformFourier(fps, DATA.DFT(ID).r*p2m, prefiltering, false);
    [~, DATA.DFT(ID).C, ~, DATA.DFT(ID).Css, DATA.DFT(ID).Cpsd]                                                     = transformFourier(fps, DATA.DFT(ID).c*p2m, prefiltering, false);
    [~, DATA.Boundary(ID).R, fss, DATA.Boundary(ID).Rss, DATA.Boundary(ID).Rpsd]                                   	= transformFourier(fps, DATA.Boundary(ID).r*p2m, prefiltering, false);
    [~, DATA.Boundary(ID).C, ~, DATA.Boundary(ID).Css, DATA.Boundary(ID).Cpsd]                                      = transformFourier(fps, DATA.Boundary(ID).c*p2m, prefiltering, false);
%     [~, DATA.imfindcirc(ID).R, fss, DATA.imfindcirc(ID).Rss, DATA.imfindcirc(ID).Rpsd]                              = transformFourier(fps, DATA.imfindcirc(ID).r*p2m, prefiltering, false);
%     [~, DATA.imfindcirc(ID).C, ~, DATA.imfindcirc(ID).Css, DATA.imfindcirc(ID).Cpsd]                                = transformFourier(fps, DATA.imfindcirc(ID).c*p2m, prefiltering, false);
    
    if DATA.settings.sectionFLAG
        for l = 1:DATA.settings.nSections
            [~, DATA.DFT(ID).section.R(:,l), fs, DATA.DFT(ID).section.Rss(:,l), DATA.DFT(ID).section.Rpsd(:,l)]   	= transformFourier(fps, DATA.DFT(ID).section.r(:,l)*p2m, prefiltering, false);
            [~, DATA.DFT(ID).section.C(:,l),  ~, DATA.DFT(ID).section.Css(:,l), DATA.DFT(ID).section.Cpsd(:,l)]   	= transformFourier(fps, DATA.DFT(ID).section.c(:,l)*p2m, prefiltering, false);
        end       
        DATA.DFT(ID).section.meanRpsd      	= mean(DATA.DFT(ID).section.Rpsd,2);
        DATA.DFT(ID).section.meanCpsd    	= mean(DATA.DFT(ID).section.Cpsd,2);
    end
    
    %% -        3. Finding standard deviation of data
%     noDrift                             = antiDrift([DATA.DFT(ID).r',DATA.DFT(ID).c']);
    DATA.Results.DFT.Rstd(ID)       	= std(DATA.DFT(ID).r*p2m);
    DATA.Results.DFT.Cstd(ID)           = std(DATA.DFT(ID).c*p2m);
    DATA.Results.Boundary.Rstd(ID)   	= std(DATA.Boundary(ID).r*p2m);
    DATA.Results.Boundary.Cstd(ID)  	= std(DATA.Boundary(ID).c*p2m);
%     DATA.Results.imfindcirc.Rstd(ID)   	= std(DATA.imfindcirc(ID).r*p2m);
%     DATA.Results.imfindcirc.Cstd(ID)  	= std(DATA.imfindcirc(ID).c*p2m);
end

%% If displacement data, then find Einstein's diffusion coefficient, D.
% D.R                                 = DATA.Results.DFT.Rstd^2/(2*W + 2*pi*W^2*);

%% * Taking mean of PSD
P                           = DATA.P(1:5:end);
for k = 1:size(P,2)
    if DATA.settings.sectionFLAG == 0
        for l = 1:5
            RpsdD(:,l)              = DATA.DFT(5*(k-1)+l).Cpsd;
            CpsdD(:,l)              = DATA.DFT(5*(k-1)+l).Cpsd;
            RpsdB(:,l)              = DATA.Boundary(5*(k-1)+l).Cpsd;
            CpsdB(:,l)              = DATA.Boundary(5*(k-1)+l).Cpsd;
        end
    else
        for l = 1:5
            RpsdD(:,l)              = DATA.DFT(5*(k-1)+l).section.meanRpsd;
            CpsdD(:,l)              = DATA.DFT(5*(k-1)+l).section.meanCpsd;
        end
    end
    DATA.meanRpsdD(:,k)         = mean(RpsdD(:,l),2);
    DATA.meanCpsdD(:,k)         = mean(CpsdD(:,l),2);
    DATA.meanRpsdB(:,k)         = mean(RpsdB(:,l),2);
    DATA.meanCpsdB(:,k)         = mean(CpsdB(:,l),2);
end

%% * Calculating force constant - Equipartition Theorem Method
RkappaEqD                	= 1e6*kB*T./(DATA.Results.DFT.Rstd).^2;
CkappaEqD                   = 1e6*kB*T./(DATA.Results.DFT.Cstd).^2;
RkappaEqB                   = 1e6*kB*T./(DATA.Results.Boundary.Rstd).^2;
CkappaEqB                   = 1e6*kB*T./(DATA.Results.Boundary.Cstd).^2;
% RkappaEqI                   = 1e6*kB*T./(DATA.Results.imfindcirc.Rstd).^2;
% CkappaEqI                   = 1e6*kB*T./(DATA.Results.imfindcirc.Cstd).^2;

fo                          = fitoptions('Method','NonlinearLeastSquares',...
                                         'Lower',       [0],...
                                         'Upper',       [10],...
                                         'StartPoint',  [0.1]);
ft                          = fittype('a*x');
c_rDEq                      = fit(DATA.P',RkappaEqD',ft,fo);
c_cDEq                      = fit(DATA.P',CkappaEqD',ft,fo);
c_rBEq                      = fit(DATA.P',RkappaEqB',ft,fo);
c_cBEq                      = fit(DATA.P',CkappaEqB',ft,fo);
% c_rIEq                      = fit(DATA.P',RkappaEqI',ft,fo);
% c_cIEq                      = fit(DATA.P',CkappaEqI',ft,fo);
    
DATA.Results.DFT.RkappaEq           = RkappaEqD;
DATA.Results.DFT.CkappaEq           = CkappaEqD;
DATA.Results.Boundary.RkappaEq      = RkappaEqB;
DATA.Results.Boundary.CkappaEq      = CkappaEqB;
% DATA.Results.imfindcirc.RkappaEq	= RkappaEqI;
% DATA.Results.imfindcirc.CkappaEq	= CkappaEqI;

for k = 1:6
    tempR                               = RkappaEqD(5*(k-1)+1:5*(k-1)+5);
    tempC                               = CkappaEqD(5*(k-1)+1:5*(k-1)+5);
    DATA.Results.DFT.RkappaEqMean(:,k)  = [mean(tempR); std(tempR)];
    DATA.Results.DFT.CkappaEqMean(:,k)  = [mean(tempC); std(tempC)];
    fcR(k)                              = 1.4*1e-6*DATA.Results.DFT.RkappaEqMean(1,1)/(2*pi*gamma)*P(k)/P(1); %median(1e-6*RkappaEqD(5*k-4:5*k)/(2*pi*gamma));
    fcC(k)                              = 1.4*1e-6*DATA.Results.DFT.CkappaEqMean(1,1)/(2*pi*gamma)*P(k)/P(1); %median(1e-6*CkappaEqD(5*k-4:5*k)/(2*pi*gamma));
end
return
%% * Calculating force constant - Power Spectrum Method
%% -        Extract PSD data and change unit from meter to nanometer
for ID = 1:30
    M                                   = 1e9;      %Size multiplikation factor (1 = m, 1e6 = um, 1e9 = nm)
    f                                   = fss(2:end);
    R(:,ID)                             = M^2*DATA.DFT(ID).Rpsd(2:end);
    C(:,ID)                             = M^2*DATA.DFT(ID).Cpsd(2:end);
%     R(R(:,ID)<1e-3,ID)                 	= 1e-3;
%     C(C(:,ID)<1e-3,ID)                	= 1e-3;
end


%% -        Determining mean PSDs
MeanR                               = mean(R,2);
MeanC                               = mean(C,2);

% Frequency range used in mean psd fitting
fitrange                            = [2 50];
ind                                 = f>=fitrange(1) & f<=fitrange(2);      %Indices corresponding to chosen fitrange

% Fitting
psdfo                               = fitoptions('Method','NonlinearLeastSquares',...
                                                 'Lower',       [M^2*D*0.1  ],...
                                                 'Upper',       [M^2*D      ],...
                                                 'StartPoint',  [M^2*D      ]);
psdft                               = fittype('D./(2.*pi^2.*(fc.^2+f.^2)) .*sinc(W.*f).^2 + E',...
                                              'independent',{'f'},'problem',{'W','E','fc'});                           
psdFitR                             = fit(f(ind),MeanR(ind),psdft,psdfo,'problem',{W,0.6,500});
psdFitC                             = fit(f(ind),MeanC(ind),psdft,psdfo,'problem',{W,0.6,800});
psdMeanR                           	= (psdFitR.D./(2.*pi^2.*(psdFitR.fc+f.^2))) .*sinc(W.*f).^2 + psdFitR.E;
psdMeanC                            = (psdFitC.D./(2.*pi^2.*(psdFitC.fc+f.^2))) .*sinc(W.*f).^2 + psdFitC.E;

% close all
% figure('position',scrnsz)
%     subplot(221),loglog(f,MeanR),hold on,loglog(f,psdMeanR,'r')
%     subplot(222),loglog(f,MeanC),hold on,loglog(f,psdMeanC,'r')
%     subplot(223),loglog(f,MeanR./psdMeanR)
%     subplot(224),loglog(f,MeanC./psdMeanC)

%% -        Removing noise peaks
% Frequency range used when finding peaks
fitrange                            = [30 Inf];
ind                                 = f>fitrange(1) & f<fitrange(2);      %Indices corresponding to chosen fitrange

% Removing Langevin baseline from spectra
MeanRsmooth                         = MeanR(ind)./psdMeanR(ind);
MeanCsmooth                         = MeanC(ind)./psdMeanC(ind);

% Finding peaks
[pksR,locsR]                        = findpeaks(MeanRsmooth,'THRESHOLD',0.04,'MINPEAKDISTANCE',100,'MINPEAKHEIGHT',2.3);
[pksC,locsC]                        = findpeaks(MeanCsmooth,'THRESHOLD',0.3,'MINPEAKDISTANCE',100,'MINPEAKHEIGHT',2.3);
fpksR                               = f(locsR) + min(fitrange);
fpksC                               = f(locsC) + min(fitrange);
wpks                                = 1;
    
indpksR                             = ones(size(f));
indpksC                             = ones(size(f));
for k = 1:size(fpksR,1)
    tempR                               = f<=fpksR(k)-wpks/2 | f>=fpksR(k)+wpks/2;
    indpksR                             = logical(indpksR.*tempR);
end
for k = 1:size(fpksC,1)
    tempC                               = f<=fpksC(k)-wpks/2 | f>=fpksC(k)+wpks/2;
    indpksC                             = logical(indpksC.*tempC);
end

close all
figure('position',scrnsz)
    subplot(221),loglog(f(ind),MeanRsmooth),hold on,loglog(fpksR,pksR,'ro')
    subplot(222),loglog(f(ind),MeanCsmooth),hold on,loglog(fpksC,pksC,'ro')
    subplot(223),loglog(f(indpksR),R(indpksR))%,hold on,loglog(f(~indpksR),R(~indpksR),'ro')
    subplot(224),loglog(f(indpksC),C(indpksC))%,hold on,loglog(f(~indpksC),C(~indpksC),'ro')


%% -        Prefit D and E
for ID = 1:30
    % Parameters to fit
    pTF                                 = logical([0 1 1]);
    
    % Frequency range used in prefitting
    fitrange                            = [fcR(ceil(ID/5)) f(end)];
    ind                                 = f>=fitrange(1) & f<=fitrange(2);      %Indices corresponding to chosen fitrange
    indR                                = logical(ind.*indpksR);
    indC                                = logical(ind.*indpksC);
    
    % Start points
    fParR                               = [0 M^2*fPar(2) M^2*fPar(3)];
    fParC                               = [0 M^2*fPar(2) M^2*fPar(3)];
    
    % Fitting model
    modelr                              = fitmodel(pTF, fParR, cPar);
    modelc                              = fitmodel(pTF, fParC, cPar);
    
    % Fitting
    prefitfo                            = fitoptions('Method','NonlinearLeastSquares',...
                                         'Lower',       [M^2*D*0.1	,     0],...
                                         'Upper',       [M^2*D*100  ,   100],...
                                         'StartPoint',  [M^2*D      , M^2*E]);
    prefitft                            = fittype('D./(2.*pi^2.*f.^2) .*sinc(W.*f).^2 + E',...
                                          'independent',{'f'},'problem',{'W'});
    prefitR                             = fit(f(indR),R(indR,ID),prefitft,prefitfo,'problem',W);
    prefitC                             = fit(f(indC),C(indC,ID),prefitft,prefitfo,'problem',W);
%     fParR(pTF)                          = expfit(modelr, favg(ind), Ravg(ind), fParR(pTF));      %Fitting exponential distributed data
%     fParC(pTF)                          = expfit(modelc, favg(ind), Cavg(ind), fParC(pTF));      %Fitting exponential distributed data
    
    
    RprefitE(ID)                        = prefitR.E;%fParR(2);
    CprefitE(ID)                        = prefitC.E;%fParC(2);
    RprefitD(ID)                        = prefitR.D;%fParR(3);
    CprefitD(ID)                        = prefitC.D;%fParC(3);
    Rprefit{ID}                      	= (RprefitD(ID)./(2.*pi^2.*f(indR).^2)) .*sinc(W.*f(indR)).^2 + RprefitE(ID);
    Cprefit{ID}                     	= (CprefitD(ID)./(2.*pi^2.*f(indC).^2)) .*sinc(W.*f(indC)).^2 + CprefitE(ID);

%     scrnsz                              = [100,100,1700,1000];%get(0,'ScreenSize');
%     figure(101);
%     set(101,'Position', scrnsz)
%     subplot(121),loglog(f(~indR), R(~indR),'bo'), hold on
%     subplot(121),loglog(f(indR), R(indR),'k'), hold on
%     subplot(121),loglog(f(indR), Rprefit{ID},'r')
%         xlim([f(1) f(end)]),xlabel('frequency, {\itf} [Hz]'),ylabel('PSD [nm^2/Hz]')
%     subplot(122),loglog(f(~indC), C(~indC),'bo'), hold on
%     subplot(122),loglog(f(indC), C(indC),'k'), hold on
%     subplot(122),loglog(f(indC), Cprefit{ID},'r')
%         xlim([f(1) f(end)]),xlabel('frequency, {\itf} [Hz]'),ylabel('PSD [nm^2/Hz]')
% 
%     figure(102)
%     subplot(121),plot(RprefitE), hold on, plot(CprefitE,'r'), title('prefit E')
%     subplot(122),plot(RprefitD), hold on, plot(CprefitD,'r'), title('prefit D')
end

meanPrefitE                         = 0.3;%mean([RprefitE(20:end),CprefitE(20:end)]);
meanPrefitD                         = mean([RprefitD(1:20),CprefitD(1:20)]);
fPar                                = [1 meanPrefitE meanPrefitD];

%% -        Fitting fc to exponential distributed PSD data
% The maximum likelihood fitting is only suitable for non-averaged data,
% using the exponential distribution characteristics of the data. Averaged
% data is gaussian distributed.
clear ind
flower                              = floor(ceil(0.2:0.2:6).^1.4); % not used
for ID = 1:30       
    % Fitting iteration
    for k = 1:2
        if k == 1
            % Frequency range used in prefitting
            fitrange                            = [9 f(end)];
            ind{ID}                             = f>=fitrange(1) & f<=fitrange(2);      %Indices corresponding to chosen fitrange   

            % Indices to fit
            indR(:,ID)                          = logical(ind{ID}.*indpksR);
            indC(:,ID)                          = logical(ind{ID}.*indpksC);
    
            % Parameters to fit
            pTF                                 = logical([1 0 1]);
            
            % Start points
            fParR                               = [fcR(ceil(ID/5)), 0, meanPrefitD];
            fParC                               = [fcC(ceil(ID/5)), 0, meanPrefitD];
            
            % Initial guess plot
            PiniR(:,ID)                         = (fParR(3)./(2.*pi^2.*(f.^2 + fParR(1)^2))) .*sinc(W.*f).^2 + fParR(2);
            PiniC(:,ID)                         = (fParC(3)./(2.*pi^2.*(f.^2 + fParC(1)^2))) .*sinc(W.*f).^2 + fParC(2);
        else
            % Frequency range used in prefitting
            fitrange                            = [9 f(end)];
            ind{ID}                             = f>=fitrange(1) & f<=fitrange(2);      %Indices corresponding to chosen fitrange   

            % Indices to fit
            indR(:,ID)                          = logical(ind{ID}.*indpksR);
            indC(:,ID)                          = logical(ind{ID}.*indpksC);
    
            % Parameters to fit
            pTF                                 = logical([1 0 1]);
            
            % Start points
            fParR                               = [fParR(1), 0, fParR(3)];
            fParC                               = [fParR(1), 0, fParC(3)];
            
            % Initial guess plot
            PiniR(:,ID)                         = (fParR(3)./(2.*pi^2.*(f.^2 + fParR(1)^2))) .*sinc(W.*f).^2 + fParR(2);
            PiniC(:,ID)                         = (fParC(3)./(2.*pi^2.*(f.^2 + fParC(1)^2))) .*sinc(W.*f).^2 + fParC(2);
        end
        
        % Fitting
        modelr                              = fitmodel(pTF, fParR, cPar);                                %Fitting model
        modelc                              = fitmodel(pTF, fParC, cPar);                                %Fitting model
        fParR(pTF)                          = expfit(modelr, f(indR(:,ID)), R(indR(:,ID),ID), fParR(pTF));           %Fitting exponential distributed data
        fParC(pTF)                          = expfit(modelc, f(indC(:,ID)), C(indC(:,ID),ID), fParC(pTF));           %Fitting exponential distributed data
        fParR                               = abs(fParR);
        fParC                               = abs(fParC);
        fitdataR(:,ID)                      = fParR;
        fitdataC(:,ID)                      = fParC;
        DATA.Results.DFT.RkappaPSD(:,ID)	= 1e6*fitdataR(1,ID)*2*pi*gamma;
        DATA.Results.DFT.CkappaPSD(:,ID)	= 1e6*fitdataC(1,ID)*2*pi*gamma;
        
        PfitR{ID}                           = (fParR(3)./(2.*pi^2.*(f(ind{ID}).^2 + fParR(1)^2))) .*sinc(W.*f(ind{ID})).^2 + fParR(2); %meanPrefitE
        PfitC{ID}                           = (fParC(3)./(2.*pi^2.*(f(ind{ID}).^2 + fParC(1)^2))) .*sinc(W.*f(ind{ID})).^2 + fParC(2); %meanPrefitE
        
        if DATA.settings.sectionFLAG
            Pfitsec                             = (fPar(3)./(2.*pi^2.*(DATA{1}.sectioned.fss.^2 + fPar(1)^2))) .*sinc(W.*DATA{1}.sectioned.fss).^2 + fPar(2);
            Sfitsec                             = sqrt(Pfitsec*DATA{1}.sectioned.fss(end)/length(Pfit));
        end
    end
end

% meanR                               = [mean(fitdataR,2), std(fitdataR,1,2)]/M^2;
% meanC                               = [mean(fitdataC,2), std(fitdataC,1,2)]/M^2;
% meanE                               = mean([meanR(1,1), meanC(1,1)]);
% meanD                               = mean([meanR(2,1), meanC(2,1)]);
% figure
% subplot(121),plot(fitdataR(1,:)/M^2),hold on,plot(fitdataC(1,:)/M^2,'r')
% subplot(122),plot(fitdataR(2,:)/M^2),hold on,plot(fitdataC(2,:)/M^2,'r')

    c_rPSD                              = fit(DATA.P',DATA.Results.DFT.RkappaPSD',ft,fo);
    c_cPSD                              = fit(DATA.P',DATA.Results.DFT.CkappaPSD',ft,fo);

for k = 1:6
    tempR                               = DATA.Results.DFT.RkappaPSD(5*k-4:5*k);
    tempC                               = DATA.Results.DFT.CkappaPSD(5*k-4:5*k);
    DATA.Results.DFT.RkappaPSDMean(:,k)	= [mean(tempR); std(tempR)];
    DATA.Results.DFT.CkappaPSDMean(:,k)	= [mean(tempC); std(tempC)];
end

% Calculating standard deviations
meanfcR = zeros(6,1); meanfcC = zeros(6,1); stdfcR = zeros(6,1); stdfcC = zeros(6,1); stdkappaR = zeros(6,1); stdkappaC = zeros(6,1);
for k = 1:6
    meanfcR(k)                      	= mean(fitdataR(1,5*k-4:5*k));
    meanfcC(k)                       	= mean(fitdataC(1,5*k-4:5*k));
    stdfcR(k)                       	= std(fitdataR(1,5*k-4:5*k));
    stdfcC(k)                           = std(fitdataC(1,5*k-4:5*k));
    stdkappaR(k)                      	= sqrt( (2*pi*gamma)^2*stdfcR(k)^2 + (2*pi*meanfcR(k))^2*sdtgamma^2 );
    stdkappaC(k)                      	= sqrt( (2*pi*gamma)^2*stdfcC(k)^2 + (2*pi*meanfcC(k))^2*sdtgamma^2 );
end

return
    %% -        PSD averaging
%     avgmethod                           = 1;
%     nbin                                = [500,500];
%     [favg, Ravg(:,ID), Rerr]         	= avgdata(nbin, avgmethod, f, R(:,ID));
%     [~   , Cavg(:,ID), Cerr]          	= avgdata(nbin, avgmethod, f, C(:,ID));
%     if ~exist('PSD','var')
%         PSD                                 = M^2*( (D./(2.*pi^2.*(fcR(ceil(ID/5))+favg.^2))) .*sinc(W.*favg).^2 + E );
%     end

%% * Plotting
maxyEq                      = max([max(RkappaEqD),max(CkappaEqD)]);
maxyPSD                     = max([max(DATA.Results.DFT.RkappaPSD),max(DATA.Results.DFT.CkappaPSD)]);
FS                          = 10;

    %% -        Standard plots
if plotFLAG
    % Plotting trace, histogram and cross-correlation
    if 0
        close all
        Mag = p2m*1e6;
        for ID = [1:5:30]%size(DATA.DFT,2)
            figure('position',scrnsz)
            subplot(3,4,2),hist(Mag*DATA.DFT(ID).r,20),title(sprintf('DFT row full, std = %1.2f nm',DATA.Results.DFT.Rstd(ID)*1e9 )),xlabel('position [um]'),xlim([-1,1]*Mag)
            subplot(3,4,3),hist(Mag*DATA.DFT(ID).c,20),title(sprintf('DFT column full, std = %1.2f nm',DATA.Results.DFT.Cstd(ID)*1e9 )),xlabel('position [um]'),xlim([-1,1]*Mag)
            subplot(3,4,6),plot(t,Mag*DATA.DFT(ID).r),xlim([t(1),t(end)]),xlabel('time [s]'),ylabel('position [um]'),ylim([-1,1]*Mag)
            subplot(3,4,7),plot(t,Mag*DATA.DFT(ID).c),xlim([t(1),t(end)]),xlabel('time [s]'),ylabel('position [um]'),ylim([-1,1]*Mag)
            subplot(3,4,1),hist(Mag*DATA.Boundary(ID).r,20),title(sprintf('Boundary row full, std = %1.2f nm',DATA.Results.Boundary.Rstd(ID)*1e9 )),xlabel('position [um]'),xlim([-1,1]*Mag)
            subplot(3,4,4),hist(Mag*DATA.Boundary(ID).c,20),title(sprintf('Boundary column full, std = %1.2f nm',DATA.Results.Boundary.Cstd(ID)*1e9 )),xlabel('position [um]'),xlim([-1,1]*Mag)
            subplot(3,4,5),plot(t,Mag*DATA.Boundary(ID).r),xlim([t(1),t(end)]),xlabel('time [s]'),ylabel('position [um]'),ylim([-1,1]*Mag)
            subplot(3,4,8),plot(t,Mag*DATA.Boundary(ID).c),xlim([t(1),t(end)]),xlabel('time [s]'),ylabel('position [um]'),ylim([-1,1]*Mag)
            
            subplot(3,4,10),plot(xcorr(Mag*DATA.DFT(ID).r,Mag*DATA.DFT(ID).r)),title('xcorr(R - DFT,Boundary)')
            subplot(3,4,11),plot(xcorr(Mag*DATA.DFT(ID).c,Mag*DATA.DFT(ID).c)),title('xcorr(C - DFT,Boundary)')
        end
        return
    end
    
    % Equipartition Method - Force constant vs. Power
    figure('Position', scrnsz)
    subplot(221),plot(DATA.P,RkappaEqD,'o','linewidth',2), hold on, plot(c_rDEq,'k'),title('row data','FontSize',FS+2), legend off
        xlim([0,1.1*max(P)]),ylim([0,1.1*maxyEq]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
        text(0.1*DATA.P(end),0.9*maxyEq,sprintf('{\\it\\kappa}_D = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_rDEq.a))
    subplot(222),plot(DATA.P,CkappaEqD,'o','linewidth',2), hold on, plot(c_cDEq,'k'), title('column data','FontSize',FS+2), legend off
        xlim([0,1.1*max(P)]),ylim([0,1.1*maxyEq]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
        text(0.1*DATA.P(end),0.9*maxyEq,sprintf('{\\it\\kappa}_D = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_cDEq.a))
    % figure
    subplot(223),plot(DATA.P,RkappaEqB,'xr','linewidth',2), hold on, plot(c_rBEq,'k'),legend off
        xlim([0,1.1*max(P)]),ylim([0,1.1*maxyEq]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
        text(0.1*DATA.P(end),0.9*maxyEq,sprintf('{\\it\\kappa}_B = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_rBEq.a))
    subplot(224),plot(DATA.P,CkappaEqB,'xr','linewidth',2), hold on, plot(c_cBEq,'k'),legend off
        xlim([0,1.1*max(P)]),ylim([0,1.1*maxyEq]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
        text(0.1*DATA.P(end),0.9*maxyEq,sprintf('{\\it\\kappa}_B = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_cBEq.a))
    % figure
%     subplot(325),plot(DATA.P,RkappaEqI,'sg','linewidth',2), hold on, plot(c_rIEq,'k'),legend off
%         xlim([0,1.1*max(P)]),ylim([0,1.1*maxy]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
%         text(0.1*DATA.P(end),0.9*maxyI,sprintf('{\\it\\kappa}_I = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_rIEq.a))
%     subplot(326),plot(DATA.P,CkappaEqI,'sg','linewidth',2), hold on, plot(c_cIEq,'k'),legend off
%         xlim([0,1.1*max(P)]),ylim([0,1.1*maxy]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
%         text(0.1*DATA.P(end),0.9*maxyI,sprintf('{\\it\\kappa}_I = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_cIEq.a))
    
	% PSD Method
    ID              = 0;
    for k = 1:6
        figure(k+1);
        set(k+1,'Name',['P = ', num2str(P(k)), 'mW'],'Position', scrnsz)
        for l = 1:5
            ID              = ID + 1;
%             R               = M^2*DATA.DFT(ID).Rpsd(2:end);
%             C               = M^2*DATA.DFT(ID).Cpsd(2:end);
            subplot(2,5,l)%, hold on
%             TTr(:,ID) = [max(DATA.DFT(ID).Rpsd) min(DATA.DFT(ID).Rpsd)];
%             TTc(:,ID) = [max(DATA.DFT(ID).Cpsd) min(DATA.DFT(ID).Cpsd)];
%         	meanDiffR(ID) = mean(DATA.DFT(ID).Cpsd./PiniR(:,ID));
%             meanDiffC(ID) = mean(DATA.DFT(ID).Cpsd./PiniC(:,ID));
%                 loglog(f, R(:,ID),'r'), hold on
                loglog(f(indR(:,ID)), smooth(R(indR(:,ID),ID),100),'k'), hold on
                loglog(f, PiniR(:,ID),'b')
                loglog(f(ind{ID}), PfitR{ID},'g')
%                 legend('data','initial guess','fit')
                xlim([f(1) f(end)]),ylim(M^2*[3e-25 2e-14]),xlabel('frequency, {\itf} [Hz]'),ylabel('PSD [nm^2/Hz]')
                
            subplot(2,5,l+5)%, hold on
%                 loglog(f, C(:,ID),'r'), hold on
                loglog(f(indR(:,ID)), smooth(C(indR(:,ID),ID),100),'k'), hold on
                loglog(f, PiniC(:,ID),'b')
                loglog(f(ind{ID}), PfitC{ID},'g')
%                 legend('data','initial guess','fit')
                xlim([f(1) f(end)]),ylim(M^2*[3e-25 2e-14]),xlabel('frequency, {\itf} [Hz]'),ylabel('PSD [nm^2/Hz]')
        end
    end
%     display(mean([meanDiffR,meanDiffC]))
%     max([TTr(1,:),TTc(1,:)])
%     min([TTr(2,:),TTc(2,:)])
    

    % PSD Method - Force constant vs. Power
    figure('Position', scrnsz)
    subplot(121)
%          plot(DATA.P,DATA.Results.DFT.RkappaPSD,'o','linewidth',0.1), hold on 
        for l = 1:size(DATA.P,2)
            text(DATA.P(l),DATA.Results.DFT.RkappaPSD(l),num2str(l))
        end
%             plot(c_rPSD,'k'), hold on,legend off
            title('row data','FontSize',FS+2)
            xlim([0,1.1*max(P)]),ylim([0,1.1*maxyPSD]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
            text(0.1*DATA.P(end),0.9*maxyPSD,sprintf('{\\it\\kappa}_D = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_rPSD.a))
    subplot(122)
        plot(DATA.P,DATA.Results.DFT.CkappaPSD,'o','linewidth',2), hold on
        plot(c_cPSD,'k'), title('column data','FontSize',FS+2), legend off
            xlim([0,1.1*max(P)]),ylim([0,1.1*maxyPSD]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
            text(0.1*DATA.P(end),0.9*maxyPSD,sprintf('{\\it\\kappa}_D = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_cPSD.a))
    
    % PSD vs. Equipartition Method - Force constant vs. Power
    figure('Position', scrnsz)
    subplot(321),errorbar(P,DATA.Results.DFT.RkappaPSDMean(1,:),DATA.Results.DFT.RkappaPSDMean(2,:),'o','linewidth',2), hold on, plot(c_rPSD,'k'),title('PSD, row data','FontSize',FS+2), legend off
        xlim([0,1.1*max(P)]),ylim([0,1.1*maxyPSD]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
        text(0.1*DATA.P(end),0.9*maxyPSD,sprintf('{\\it\\kappa}_D = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_rPSD.a))
    subplot(322),errorbar(P,DATA.Results.DFT.CkappaPSDMean(1,:),DATA.Results.DFT.CkappaPSDMean(2,:),'o','linewidth',2), hold on, plot(c_cPSD,'k'), title('PSD, column data','FontSize',FS+2), legend off
        xlim([0,1.1*max(P)]),ylim([0,1.1*maxyPSD]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
        text(0.1*DATA.P(end),0.9*maxyPSD,sprintf('{\\it\\kappa}_D = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_cPSD.a))
    subplot(323),errorbar(P,DATA.Results.DFT.RkappaEqMean(1,:),DATA.Results.DFT.RkappaEqMean(2,:),'o','linewidth',2), hold on, plot(c_rDEq,'k'),title('Equipartition, row data','FontSize',FS+2), legend off
        xlim([0,1.1*max(P)]),ylim([0,1.1*maxyPSD]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
        text(0.1*DATA.P(end),0.9*maxyPSD,sprintf('{\\it\\kappa}_D = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_rDEq.a))
    subplot(324),errorbar(P,DATA.Results.DFT.CkappaEqMean(1,:),DATA.Results.DFT.CkappaEqMean(2,:),'o','linewidth',2), hold on, plot(c_cDEq,'k'), title('Equipartition, column data','FontSize',FS+2), legend off
        xlim([0,1.1*max(P)]),ylim([0,1.1*maxyPSD]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
        text(0.1*DATA.P(end),0.9*maxyPSD,sprintf('{\\it\\kappa}_D = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_cDEq.a))
    subplot(325),plot(fitdataR(1,:),'r'), hold on,plot(fitdataC(1,:)),legend('row','coloumn')
        xlabel('measurement number, [#]','FontSize',FS),ylabel('corner frequency, {\it\kappa} [Hz]','FontSize',FS)
    subplot(326),plot(fitdataR(3,:),'r'), hold on,plot(fitdataC(3,:)),legend('row','coloumn')
        xlabel('measurement number, [#]','FontSize',FS),ylabel('diffusion constant, {\itD} [m^2/s]','FontSize',FS)
    
    
    
    %% -        Plotting figures for the paper
    % Kappa vs power
    FS = 10;
    LW = 1;
    MS = 3;
    c_cPSD3um       = c_cPSD;
    c_rPSD3um       = c_rPSD;
    stdkappaR3um    = stdkappaR;
    stdkappaC3um    = stdkappaC;
    KAPPAr3um       = DATA.Results.DFT.RkappaPSDMean(1,:);
    KAPPAc3um       = DATA.Results.DFT.CkappaPSDMean(1,:);
    left = 5; bottom = 5; width = 21; height = 10;
    fig = figure(101);
    set(fig,'units','centimeters ','Position', [left, bottom, width, height])
    subplot(121),errorbar(P,KAPPAc3um,1e6*stdkappaC,'o','linewidth',LW,'MarkerSize',MS), hold on, plot([0,1.1*max(P)], [0,1.1*max(P)]*c_cPSD.a,'k','linewidth',LW), title('x-axis','FontSize',FS), legend off
        xlim([0,1.1*max(P)]),ylim([0,150]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
%         text(0.1*DATA.P(end),120,sprintf('{\\it\\kappa}_D = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_cPSD.a))
    subplot(122),errorbar(P,KAPPAr3um,1e6*stdkappaR,'o','linewidth',LW,'MarkerSize',MS), hold on, plot([0,1.1*max(P)], [0,1.1*max(P)]*c_rPSD.a,'k','linewidth',LW),title('y-axis','FontSize',FS), legend off
        xlim([0,1.1*max(P)]),ylim([0,150]),xlabel('power, {\itP} [mW]','FontSize',FS),ylabel('force constant, {\it\kappa} [pN/um]','FontSize',FS)
%         text(0.1*DATA.P(end),120,sprintf('{\\it\\kappa}_D = {\\ita}*{\\itP}, {\\ita} = %1.3f pN/(um*mW)',c_rPSD.a))
    
%% Einstein Diffusions Constant results
Dmeasured = ['D measured: ',num2str(mean([fitdataR(3,:),fitdataC(3,:)])/M^2)];
Dcalculated = ['D calculated: ',num2str(Dcalc)];
display(Dmeasured)
display(Dcalculated)

%% Saving

% save('PaperData3um','KAPPAr3um','KAPPAc3um','stdkappaR3um','stdkappaC3um','c_rPSD3um','c_cPSD3um')

return
%% PSD
    % PSD
    fig = figure(102);
    ID = 13;
    set(fig,'units','centimeters ','Position', [left, bottom, width, height])
    subplot(121)
        loglog(f, R(:,ID),'k'), hold on
        loglog(f(indR(:,ID)), R(indR(:,ID),ID),'b'), hold on
            xlim([f(1) f(end)]),ylim([10^(-5) 10^3]),xlabel('frequency [Hz]'),ylabel('distance squared per frequency [a. u.]')
%             title('3um bead trapped at 44mW')
    subplot(122)
    find = f(indR(:,ID));
        loglog(f(indR(:,ID)), smooth(R(indR(:,ID),ID),100),'b'), hold on
        loglog(f(ind{ID}), PfitR{ID},'r')
            xlim([find(1) f(end)]),ylim([10^(-5) 10^3]),xlabel('frequency [Hz]'),ylabel('distance squared per frequency [a. u.]')
%             title('3um bead trapped at 44mW')
    
    %% -        Plotting figures for group talk
    if 0
        close all
        FS = 14;
        Mag = 1;%p2m*1e6;
        for ID = 21%size(DATA.DFT,2)
            figure('position',scrnsz)
            subplot(2,2,4),hist(Mag*DATA.DFT(ID).r,20),title(sprintf('DFT, std = %1.2f nm',DATA.Results.DFT.Rstd(ID)*1e9),'FontSize',FS)
                xlabel('y-axis [pixels]','FontSize',FS),xlim([-1,1]*Mag),ylabel('counts [#]','FontSize',FS)
            subplot(2,2,3),hist(Mag*DATA.DFT(ID).c,20),title(sprintf('DFT, std = %1.2f nm',DATA.Results.DFT.Cstd(ID)*1e9),'FontSize',FS)
                xlabel('x-axis [pixels]','FontSize',FS),xlim([-1,1]*Mag),ylabel('counts [#]','FontSize',FS)
            subplot(2,2,2),hist(Mag*DATA.Boundary(ID).r,20),title(sprintf('Boundary, std = %1.2f nm',DATA.Results.Boundary.Rstd(ID)*1e9),'FontSize',FS)
                xlabel('y-axis [pixels]','FontSize',FS),xlim([-1,1]*Mag),ylabel('counts [#]','FontSize',FS)
            subplot(2,2,1),hist(Mag*DATA.Boundary(ID).c,20),title(sprintf('Boundary, std = %1.2f nm',DATA.Results.Boundary.Cstd(ID)*1e9),'FontSize',FS)
                xlabel('x-axis [pixels]','FontSize',FS),xlim([-1,1]*Mag),ylabel('counts [#]','FontSize',FS)
        end
    end
    
    return
    %% -        Plotting figures for Kishan
% Trace and Hist
close all
Mag     = p2m*1e6;
k       = 0;
for ID = [2:5:30]%size(DATA.DFT,2)
    k       = k + 1;
%     figure('position',scrnsz)
    subplot(6,2,2*k-1),hist(Mag*DATA.DFT(ID).r,20),title(sprintf('DFT row full, std = %1.2f nm',DATA.Results.DFT.Rstd(ID)*1e9 )),xlabel('position [um]'),xlim([-1,1]*Mag)
    subplot(6,2,2*k),hist(Mag*DATA.DFT(ID).c,20),title(sprintf('DFT column full, std = %1.2f nm',DATA.Results.DFT.Cstd(ID)*1e9 )),xlabel('position [um]'),xlim([-1,1]*Mag)
%     subplot(3,2,3),plot(t,Mag*DATA.DFT(ID).r),xlim([t(1),t(end)]),xlabel('time [s]'),ylabel('position [um]'),ylim([-1,1]*Mag)
%     subplot(3,2,4),plot(t,Mag*DATA.DFT(ID).c),xlim([t(1),t(end)]),xlabel('time [s]'),ylabel('position [um]'),ylim([-1,1]*Mag)
%     subplot(3,2,5),plot(xcorr(Mag*DATA.DFT(ID).r,Mag*DATA.DFT(ID).r)),title('xcorr(R - DFT,Boundary)')
%     subplot(3,2,6),plot(xcorr(Mag*DATA.DFT(ID).c,Mag*DATA.DFT(ID).c)),title('xcorr(C - DFT,Boundary)')
end

% PSD
ID              = 3;
figure(102);
set(102,'Position', scrnsz)
for k = 1:6
%     subplot(2,3,k)%, hold on
%     loglog(f, M^2*DATA.DFT(ID).Rpsd,'r'), hold on
%     loglog(f(indR(:,ID)), M^2*DATA.DFT(ID).Rpsd(indR(:,ID)),'k'), hold on
%     loglog(f, PiniR(:,ID),'b')
%     loglog(f(ind), PfitR(:,ID),'g')
%         xlim([f(1) f(end)]),ylim(M^2*[3e-25 2e-14])
%         xlabel('frequency, {\itf} [Hz]'),ylabel('PSD [nm^2/Hz]'),title(sprintf('P = %1.0f mW',P(k)))
    
    subplot(2,3,k)%, hold on
    loglog(f, M^2*DATA.DFT(ID).Cpsd,'r'), hold on
    loglog(f(indR(:,ID)), M^2*DATA.DFT(ID).Cpsd(indR(:,ID)),'k'), hold on
    loglog(f, PiniC(:,ID),'b')
    loglog(f(ind), PfitC(:,ID),'g')
        xlim([f(1) f(end)]),ylim(M^2*[3e-25 2e-14])
        xlabel('frequency, {\itf} [Hz]'),ylabel('PSD [nm^2/Hz]'),title(sprintf('P = %1.0f mW',P(k)))
    
    ID              = ID + 5;
end

% PSD averaging
figure(103);
set(103,'Position', scrnsz)
for ID = 1:15;
    avgmethod                   = 2;
    nbin                        = [1000,200];
    ffull                       = f(2:end);
    xpfull                      = M^2*DATA.DFT(ID).Rpsd(2:end);
    [favg, xavg, xerr]          = avgdata(nbin, avgmethod, ffull, xpfull);

    subplot(3,5,ID)
    loglog(favg, xavg,'ko'), hold on
    loglog(f(ind), PfitR(:,ID),'r')
        xlim([f(1) f(end)]),xlabel('frequency, {\itf} [Hz]'),ylabel('PSD [nm^2/Hz]')
end 

    %% -        Extra plots
    return
    % Comparison of DFT and Boundary PSDs
    for k = 1:6
        figure(k+9)
        set(k+9,'Position', scrnsz)
        subplot(3,1,1),loglog(f, DATA.meanRpsdD(:,k)),title('DFT psd')
        subplot(3,1,2),loglog(f, DATA.meanRpsdB(:,k)),title('Boundary psd')
        subplot(3,1,3),loglog(f, DATA.meanRpsdD(:,k)./DATA.meanRpsdB(:,k)),hold on
        subplot(3,1,3),loglog(f,ones(size(f)),'r'),title('DFT psd/Boundary psd')
    end
    
    
%     figure(2)
%     subplot(1,2,1)%, hold on
%         loglog(f,DATA.DFT(ID).Rpsd,'k')
%     subplot(1,2,2)%, hold on
%         loglog(f,DATA.DFT(ID).Cpsd,'k')
%     subplot(1,2,1), hold on
%         loglog(f1(ind),Pfit,'r')
%         xlim([1e0 1e7])
%         ylim([1e-8 1e8])
%         set(gca,'YTick',[1e-6, 1e-3, 1e0, 1e3, 1e6])
%         title('Full dataset')
%         ylabel('Power spectral density')
%     subplot(3,2,2)
%         semilogy(f1,Sss1,'k')
%     subplot(3,2,2), hold on
%         semilogy(f1(ind),Sfit,'r')
%         xlim([1e4,3e6])
%         ylim([1e-5 1e2])
%         set(gca,'YTick',[1e-4, 1e-2, 1e0, 1e2])
%         ylabel('Fourier Transform')
%     subplot(3,2,3), hold on
%         plot(f1(ind),Sss1(ind)./Sfit,'k')
%         ylabel('Fourier Transform / fit')
    
    if DATA.settings.sectionFLAG ==3
        figure(3)
        subplot(3,1,1)%, hold on
            loglog(DATA{1}.sectioned.fss,PSD_smooth,'k')
        subplot(3,1,1), hold on
            loglog(DATA{1}.sectioned.fss,Pfitsec,'r')
            xlim([1e0 1e7])
            ylim([1e-8 1e8])
            set(gca,'YTick',[1e-6, 1e-3, 1e0, 1e3, 1e6])
            title(sprintf('Sectioned dataset, n = %2.0u',n))
            ylabel('Power spectral density')
        subplot(3,1,2)
            semilogy(DATA{1}.sectioned.fss,Sss_smooth,'k')  
        subplot(3,1,2), hold on
            semilogy(DATA{1}.sectioned.fss,Sfitsec,'r')
            xlim([1e4,3e6])
            ylim([1e-5 1e2])
            set(gca,'YTick',[1e-4, 1e-2, 1e0, 1e2])
            ylabel('Fourier Transform')
        subplot(3,1,3), hold on
            plot(DATA{1}.sectioned.fss,Sss_smooth./Sfitsec,'r')
            ylabel('Fourier Transform / fit')
    end
end

