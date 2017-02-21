close all
clc

if 1 %Switches between tracking and analysis
%%
for fileID = 2
clearvars -except fileID
dataFolder                 	= 'Jonny and Claire - Data\2014_05_22\800nmGRINlens';
FOLDER                      = [cd, filesep, dataFolder, filesep, num2str(fileID+1),'um'];     	%Full path to folder containing the data files
TYPE                        = 'mat';                                   	%File format
SIZE                        = [30, 1000]*10^6;                        	%Range of files size in MB
[FILES]                   	= filenames(FOLDER,TYPE,SIZE,false);

%%
for ID = 1:30
if 1
    clearvars -except DATA ID FOLDER FILES settings fileID F nFrames
    filename                    = FILES{ID,:};
    DATA.P(ID)                  = str2double(filename(5:6));
    load([FOLDER,filesep,filename,'.mat'])
    if ID == 1
        p2m                         = 0.1233*1e-6;              	%Pixel to meter
        T                           = 293.15;                      	%Temperature [K]
        F                           = 1;
        nFrames                     = settings.nFrames;
        fps                         = 400;%settings.fps;
        settings.p2m                = p2m;
        settings.T                  = T;
        settings.F                	= F;
        settings.fps                = fps;
        settings.d                  = 4*1e-6;
        DATA.settings               = settings;
    end
    
    %% Show movie
%     figure
%     for k = 1:40000
%         frame = mov(:,:,k);
%         imagesc(frame),drawnow
%         subplot(221),imagesc(frame(16:105,7:96)),axis image,colormap(gray),ylabel('11 um beads','Fontsize',14),title('Original images','Fontsize',14)
%         subplot(222),imagesc(mov(16:105,7:96,k)),axis image,colormap(gray),drawnow,title('Filtered images','Fontsize',14)
%         subplot(223),imagesc(frame(11+10:100+10,7+15:96+15)),axis image,colormap(gray),ylabel('5 um beads','Fontsize',14)
%         subplot(224),imagesc(mov(11+10:100+10,7+15:96+15,k)),axis image,colormap(gray),drawnow
%     end
end

noCircleFlag                = 0; %counters for possible errors
multipleCircleFlag          = 0; %counters for possible errors
displacementFlag            = 1;
[Xin,Yin]                   = meshgrid(1:size(mov,1));
[Xout,Yout]                 = meshgrid(1:1/F:size(mov,1));
if displacementFlag
    nmax                        = nFrames-1;
else
    nmax                        = nFrames;
end

%Preallocation
DATA.DFT(ID).error          = zeros(nmax,1);
DFT.r                       = zeros(nmax,1);
DFT.c                       = zeros(nmax,1);
loadtime                    = zeros(nmax,1);
dfttime                     = zeros(nmax,1);
boundarytime                = zeros(nmax,1);
imfindcirctime           	= zeros(nmax,1);

for n = 1:nmax
    %load frame
    tic
    if displacementFlag
        frame0                      = mov(:,:,n);
        frame1                      = mov(:,:,n+1);
    else
        frame0                      = mov(:,:,n);
        if F ~= 1
            frame1                      = interp2(Xin,Yin,frame0,Xout,Yout,'linear');
        else
            frame1                      = frame0;
        end
    end
    loadtime(n)                	= toc;
    
    %% DFT Tracking
    tic
    if 1
        % Preparing images
        if displacementFlag
            imR                         = frame0 - mean(frame0(:));
            imT                         = frame1 - mean(frame1(:));
        else
            if n == 1
                imR                         = frame1 - mean(frame1(:));
            end
            imT                         = frame1 - mean(frame1(:));
        end
        %                 imagesc(imT),colormap(gray),axis image xy, drawnow
        
        % Determining translation
        usfac                       = 200; %Resolution is: [1/usfac] pixels
        [output, ~]                 = dftregistration_v1(fft2(imR),fft2(imT),usfac);
        DATA.DFT(ID).error(n)       = output(1);
        DFT.r(n)                    = output(3);
        DFT.c(n)                    = output(4);
    end
    dfttime(n)                	= toc;
    
    %% Boundary Tracking
    tic
    if 1
        I                       = frame1;
        
        % Color to black and white
        mu                      = mean(I(:));
        sd                      = std(I(:));
        th                      = 0.10;%0.994;%mu+0.6*sd;
        bw0                     = im2bw(I,th);
        
        % Cleaning up image
        n1                      = 0;
        bw1                     = imclose(bw0,strel('disk', n1));
        n2                      = 0;
        bw2                     = bwareaopen(bw1,n2);
        bw3                     = imfill(bw2,'holes');
        
        % Get the boundary
        [B1]                    = bwboundaries(bw3,'noholes');
        Br{n}                   = B1{1}(:,1);
        Bc{n}                   = B1{1}(:,2);
        
        % Plotting
        if 1
            subplot(231),imagesc(I),axis image
            subplot(232),imagesc(bw0),axis image
            subplot(233),imagesc(bw1),axis image
            subplot(234),imagesc(bw2),axis image
            subplot(235), imagesc(bw3),axis image, colormap(gray)
            hold on
            plot(Bc{n},Br{n},'r')
            plot(mean(Bc{n}),mean(Br{n}),'rx','MarkerSize',10)
            drawnow
        end
        
        % Figure plot
        if 0
            imagesc(I),axis image
            set(gca,'xtick',[],'ytick',[]),colormap('gray')
            hold on
            plot(Bc{n},Br{n},'r')
            plot(mean(Bc{n}),mean(Br{n}),'rx','MarkerSize',10)
        end
    end
    boundarytime(n)        	= toc;
    return
    %% Using imfindcircles function
    tic
    if 0
        frame                       = abs((1-frame1).^1); % F=1(^1.1),  F=2(.^0.5)  ^0.9  ^0.9
        [centers, radii, metric]    = imfindcircles(frame,[30 45],'ObjectPolarity','bright'); % F=1[10 20],  F=2[20 30],  [47 70]
        
        %check for anomalies
        if isempty(centers) == 1
            figure();imagesc(frame);
            title(strcat('frame # ',num2str(n),': no circles found'));axis image;colormap gray;
            drawnow;
            pause();
            close gcf
            CoM(n,:)                    = [NaN,NaN];
            R(n,:)                      = NaN;
            GoT(n,:)                    = NaN;
            noCircleFlag                = noCircleFlag+1;
        elseif max(size(centers) ~= [1,2])
            figure();imagesc(frame);
            title(strcat('frame # ',num2str(n),': more than one circle found, choosing green circle'));axis image;colormap gray;
            viscircles(centers, radii,'EdgeColor','r');
            viscircles(centers(1,:), radii(1),'EdgeColor','g');
            drawnow;
            pause();
            close gcf
            CoM(n,:)                    = centers(1,:);
            R(n,:)                      = radii(1,:);
            GoT(n,:)                    = metric(1,:);
            multipleCircleFlag          = multipleCircleFlag+1;
        else
%             imagesc(frame);
%             title(strcat('frame # ',num2str(n),': found one circle'));axis image;colormap gray;
%             viscircles(centers(1, :), radii(1),'EdgeColor','g');
%             drawnow;
            
            CoM(n,:)                    = centers;
            R(n,:)                      = radii;
            GoT(n,:)                    = metric;
        end
    end
    imfindcirctime(n)         	= toc;
    
    % Progress
    if n/round(nFrames/100) == round(n/round(nFrames/100))
        clc
        display(['fileID = ',num2str(fileID),', ID = ',num2str(ID),', frame # = ',num2str(n),', ',num2str(100*n/nFrames),'%'])
        display(['Time per cycle: Loading frame = ',num2str(mean(loadtime(1:n))),'s, DFT = ',num2str(mean(dfttime(1:n))),'s, Boundary = ',num2str(mean(boundarytime(1:n))),'s, imfindcircles = ',num2str(mean(imfindcirctime(1:n))),'s'])
    end
end

%% Restructure data
if 1
    for l = 1:nmax
        c(l)                  	= mean(Bc{l});
        r(l)                	= mean(Br{l});
    end
    tempr                   = (r - mean(r))./F;
    tempc                   = (c - mean(c))./F;
    if displacementFlag
        DATA.Boundary(ID).r     = diff(tempr);
        DATA.Boundary(ID).c     = diff(tempc);
    else
        DATA.Boundary(ID).r     = tempr;
        DATA.Boundary(ID).c     = tempc;
    end
    
    DATA.DFT(ID).r          = (DFT.r - mean(DFT.r))./F;
    DATA.DFT(ID).c          = (DFT.c - mean(DFT.c))./F;
    
%     DATA.imfindcirc(ID).r  	= (CoM(:,2)-mean(CoM(:,2)))./F;
%     DATA.imfindcirc(ID).c  	= (CoM(:,1)-mean(CoM(:,1)))./F;
end

end
    save([FOLDER,filesep,'DATA_F1_',num2str(fileID+1),'um_displacement.mat'], 'DATA')
end
return
end


%% Analysis
t                    	= 0:1/DATA.settings.fps:(DATA.settings.nFrames-1)/DATA.settings.fps;

for ID = 1:30
    %% Equipartition
    DFT                     = [DATA.DFT(ID).r           ;  	DATA.DFT(ID).c          ];
    Boundary                = [DATA.Boundary(ID).r      ;	DATA.Boundary(ID).c     ];
    imfindcirc              = [DATA.imfindcirc(ID).r'   ;   DATA.imfindcirc(ID).c'  ];
    DFTstd.full(ID,:)    	= std(DFT');
    
    %% Power spectra
%     PSDr = zeros(418,1);
%     PSDc = PSDr;
%     for k = 1:3
%         [U,S,Uss,Sss,PSD] = transformFourier(t(1:3:end-2),DFT(1,k:3:end-(3-k)),0,0);
%         PSDr = PSDr + PSD;
%         [U,S,Uss,Sss,PSD] = transformFourier(t(1:3:end-2),DFT(2,k:3:end-(3-k)),0,0);
%         PSDc = PSDc + PSD;
%     end
%     subplot(121),loglog(Uss,PSDr)
%     subplot(122),loglog(Uss,PSDc)

    %% Plotting
%     figure
%     subplot(2,2,1),hist(DFT(1,:),20),title(sprintf('DFT row full, std = %1.2f pixels',DFTstd.full(ID,1))),xlabel('position [pixels]'),xlim([-1,1])
%     subplot(2,2,2),hist(DFT(2,:),20),title(sprintf('DFT column full, std = %1.2f pixels',DFTstd.full(ID,2))),xlabel('position [pixels]'),xlim([-1,1])
%     subplot(2,2,3),plot(t,DFT(1,:)),xlim([t(1),t(end)]),xlabel('time [s]'),ylabel('position [pixels]'),ylim([-1,1])
%     subplot(2,2,4),plot(t,DFT(2,:)),xlim([t(1),t(end)]),xlabel('time [s]'),ylabel('position [pixels]'),ylim([-1,1])
%     
%     figure
% 	id = 2;
%     int = [-1.5,1.5];
%     subplot(3,3,1),hist(DFT(id,:),20),title(sprintf('DFT, std = %1.2f',std(DFT(id,:)))),xlim(int)
%     subplot(3,3,2),hist(Boundary(id,:),20),title(sprintf('Boundary, std = %1.2f',std(Boundary(id,:)))),xlim(int)
%     subplot(3,3,3),hist(imfindcirc(id,:),20),title(sprintf('imfindcirc, std = %1.2f',std(imfindcirc(id,:)))),xlim(int)
%     subplot(3,3,4),plot(t,DFT(id,:)),xlim([t(1),t(end)]),ylim(int)
%     subplot(3,3,5),plot(t,Boundary(id,:)),xlim([t(1),t(end)]),ylim(int)
%     subplot(3,3,6),plot(t,imfindcirc(id,:)),xlim([t(1),t(end)]),ylim(int)
%     subplot(3,3,7),plot(xcorr(DFT(id,:),Boundary(id,:))),title('xcorr(DFT,Boundary)')
%     subplot(3,3,8),plot(xcorr(Boundary(id,:),imfindcirc(id,:))),title('xcorr(Boundary,imfindcirc)')
%     subplot(3,3,9),plot(xcorr(DFT(id,:),imfindcirc(id,:))),title('xcorr(DFT,imfindcirc)')
end

%% Calculating force constant
kappa_r                 = 1e6*kB*T./(p2m*DFTstd.full(:,1)).^2;
kappa_c                 = 1e6*kB*T./(p2m*DFTstd.full(:,2)).^2;

ft                      = fittype('a*x');
c_r                     = fit(P',kappa_r,ft);
c_c                     = fit(P',kappa_c,ft);

maxy                    = max([max(kappa_r),max(kappa_c),max(kappa_r),max(kappa_c)]);
figure
subplot(121),plot(P,kappa_r,'o','linewidth',2), hold on, plot(c_r),title('\kappa_{row} full data'),xlim([0,1.1*max(P)]),ylim([0,1.1*maxy]),xlabel('power [mW]'),ylabel('force constant [pN/um]'), legend off
    text(0.1*P(end),0.9*maxy,sprintf('y = a*x, a = %1.3f',c_r.a))
subplot(122),plot(P,kappa_c,'o','linewidth',2), hold on, plot(c_c), title('\kappa_{column} full data'),xlim([0,1.1*max(P)]),ylim([0,1.1*maxy]),xlabel('power [mW]'),ylabel('force constant [pN/um]'), legend off
    text(0.1*P(end),0.9*maxy,sprintf('y = a*x, a = %1.3f',c_c.a))
    
return
%% Plotting
figure
subplot(221),plot(t,c),ylim([-1,1])
subplot(222),plot(t,r),ylim([-1,1])
subplot(223),plot(t,F.c),ylim([-1,1])
subplot(224),plot(t,F.r),ylim([-1,1])

figure
subplot(221),hist(B.c,20),xlim([-1,1])
subplot(222),hist(B.r,20),xlim([-1,1])
subplot(223),hist(F.c,20),xlim([-1,1])
subplot(224),hist(F.r,20),xlim([-1,1])

figure
subplot(121),plot(xcorr(F.c,c))
subplot(122),plot(xcorr(F.r,r))

figure
[U,S,Uss,Sss,PSD] = transformFourier(t,r,0,0);
subplot(121),loglog(Uss,PSD),hold on
[U,S,Uss,Sss,PSD] = transformFourier(t,F.r,0,0);
subplot(121),loglog(Uss,PSD,'r')
[U,S,Uss,Sss,PSD] = transformFourier(t,c,0,0);
subplot(122),loglog(Uss,PSD),hold on
[U,S,Uss,Sss,PSD] = transformFourier(t,F.c,0,0);
subplot(122),loglog(Uss,PSD,'r')

%% Analysis of beads at ~400Hz
return
x(:,3) = CoM(:,1);
y(:,3) = CoM(:,2);
return
PSDx = zeros(1501,1);
PSDy = PSDx;
for k = 1:3
    [U,S,Uss,Sss,PSD] = transformFourier(t,x(:,k),0,0);
    PSDx = PSDx + PSD;
    [U,S,Uss,Sss,PSD] = transformFourier(t,y(:,k),0,0);
    PSDy = PSDy + PSD;
end
subplot(121),loglog(Uss,PSD)
subplot(122),loglog(Uss,PSD)

figure,plot(t,CoM(:,2))
