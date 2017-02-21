close all
clear
clc

imSize              = 50; % 11um: 80, 5um: 60
scrnsz            	= [100,100,1700,1000];%get(0,'ScreenSize');
%% Making a list of individual file names in the "Data" folder
cf                  = cd;
temp                = regexp(cf,filesep);
BaseFolder          = cf(1:temp(end)-1);
DataFolder          = 'Test_Data';             %Folder containing the data files
FOLDER              = ['C:\Users\marv\Documents\Work\130601-150531 Postdoc Villum Fonden\130601-141130 Postdoc St Andrews\Open Projects\Jonny&Claire - SIM  fiber trap\Jonny and Claire - Data\2014_05_22\800nmGRINlens\4um'];%, filesep, DataFolder];              	%Full path to folder containing the data files
TYPE                = 'avi';                                   	%File format
SIZE                = [1, 3000]*10^6;                        	%Range of files size in MB
FILES            	= filenames(FOLDER,TYPE,SIZE,false);

PPall               = zeros(20,30);
for movnr = 1:size(FILES,1) %17 used in fiber trap paper
    %% Loading movie
    inputFile           = [FOLDER,'\',FILES{movnr,:},'.avi'];
    outputFile          = [FOLDER,'\',FILES{movnr,:},'.mat'];
    inputVideo          = VideoReader(inputFile);
    nFrames             = inputVideo.NumberOfFrames;
    rFrames             = 2:nFrames;
    settings.nFrames    = size(rFrames,2);
    settings.fps       	= inputVideo.FrameRate;
    T                   = (settings.nFrames-1)/settings.fps;
    f                   = linspace(1/T,settings.fps,settings.nFrames);
    
    %% Preallocation
    clear mov3 movBG0 movBG
    mov0                        = zeros(inputVideo.Height,inputVideo.Width,settings.nFrames);
    mov1                        = mov0;
%     mov2                        = zeros(imSize,imSize,settings.nFrames);
%     mov4                        = mov2;
    
    %% Creating movie array
    for k = rFrames
        kk                  = 1 + k-rFrames(1);
        frame           	= read(inputVideo, k);
        frame            	= double(frame)/255;%double(rgb2gray(frame))/255;
        mov0(:,:,kk)    	= frame;
%             imagesc(frame), axis image
%             drawnow,pause(0.03)
        if k/10 == round(k/10)
            clc
            display(['Movie nr.:', num2str(movnr), sprintf(', Progress: %1.1f%%', 100*kk/settings.nFrames)])
        end
    end
    
    
     %% Removing Slope background
    if 0
        me                  = mean(mov0,3);
        se                  = std(mov0,1,3);
        [X,Y]               = meshgrid(1:size(me,2),1:size(me,1));
        x                   = reshape(X,numel(me),1);
        y                   = reshape(Y,numel(me),1);
        z                   = reshape(me,numel(me),1);
        
        fo                  = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',       [  -Inf, -Inf,    -Inf],...
            'Upper',       [Inf,  Inf,  Inf],...
            'StartPoint',  [ 0  0     0    ]);
        ft                  = fittype('a*x+b*y+offset',...
            'dependent',{'z'},'independent',{'x','y'});
        c0                  = fit([x,y],z,ft,fo);
        zfit0               = c0.a*x+c0.b*y+c0.offset';
%         c1                  = fit([x(zfit0-z>0.0),y(zfit0-z>0.0)],z(zfit0-z>0.0),ft,fo);
%         zfit1               = c1.a*x+c1.b*y+c1.offset';
%         c                   = fit([x(zfit1-z>0.0),y(zfit1-z>0.0)],z(zfit1-z>0.0),ft,fo);
%         zfit                = c.a*x+c.b*y+c.offset';
        zfit0               = reshape(zfit0,size(me,1),size(me,2));
%         zfit1               = reshape(zfit1,size(me,1),size(me,2));
%         zfit                = reshape(zfit,size(me,1),size(me,2));
%             figure
%                 subplot(131),surf(zfit0,'EdgeColor','none')
%                 hold on,surf(me,'EdgeColor','none')
%                 subplot(132),surf(zfit1,'EdgeColor','none')
%                 hold on,surf(me,'EdgeColor','none')
%                 subplot(133),surf(zfit,'EdgeColor','none')
%                 hold on,surf(me,'EdgeColor','none')
        for k = 1:settings.nFrames
            mov1(:,:,k)              	= mov0(:,:,k) - zfit0;
%                 subplot(1,2,1),imagesc(mov0(:,:,k)), axis image,subplot(1,2,2),imagesc(mov1(:,:,k) ), axis image
%                 drawnow,pause(0.03)
        end
    else
        mov1                      	= mov0;
    end
    
    %% Removing Gaussian background
    if 0
        me                  = mean(mov1,3);
        se                  = std(mov1,1,3);
        [X,Y]               = meshgrid(1:size(me,2),1:size(me,1));
        x                   = reshape(X,numel(me),1);
        y                   = reshape(Y,numel(me),1);
        z                   = reshape(me,numel(me),1);
        
        fo                  = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',       [  0, -100,    10,   -Inf,  -Inf],...
            'Upper',       [100,  100,  10000,   Inf,   Inf],...
            'StartPoint',  [0.4     0      50     30     30]);
        ft                  = fittype('A*exp(-((x-x0)^2+(y-y0)^2)/(2*sigma^2))+offset',...
            'dependent',{'z'},'independent',{'x','y'});
        c0                  = fit([x,y],z,ft,fo);
        zfit0               = c0.A*exp(-((x-c0.x0).^2+(y-c0.y0).^2)./(2*c0.sigma^2))+c0.offset;
        c1                  = fit([x(zfit0-z<0.05),y(zfit0-z<0.05)],z(zfit0-z<0.05),ft,fo);
        zfit1               = c1.A*exp(-((x-c1.x0).^2+(y-c1.y0).^2)./(2*c1.sigma^2))+c1.offset;
        c                   = fit([x(zfit1-z<0.05),y(zfit1-z<0.05)],z(zfit1-z<0.05),ft,fo);
        zfit                = c.A*exp(-((x-c.x0).^2+(y-c.y0).^2)./(2*c.sigma^2))+c.offset;
        zfit                = reshape(zfit,size(me,1),size(me,2));
            figure,surf(zfit,'EdgeColor','none')
            hold on,surf(me,'EdgeColor','none')
        for k = 1:settings.nFrames
            mov1(:,:,k)              	= mov1(:,:,k) ./ zfit;
%                 subplot(1,2,1),imagesc(mov0(:,:,k)), axis image,subplot(1,2,2),imagesc(mov1(:,:,k) ), axis image
%                 drawnow,pause(0.03)
        end
    elseif sum(mov1(:)) == 0
        mov1                      	= mov0;
    end
    
    %% Cropping
    if 1
        if DataFolder(end) == '5'
            me1                         = abs(mean(mov1,3).^0.5);
            B                           = 20;
            [center, radius, metric]    = imfindcircles(abs(me1(1+B:end-B,1+B:end-B)),[10 20],'ObjectPolarity','bright');
        elseif DataFolder(end) == '11'
            me1                         = mean(mov1,3);
            B                           = 0;
            [center, radius, metric]    = imfindcircles(me1,[25 35],'ObjectPolarity','dark');
        else
            me1                         = abs(mean(mov1,3).^0.8); % 4um: abs(mean(mov1,3).^0.8);
            me2                         = me1;
            me2(me2>0.0007)          	= 0.0007;
            B                           = 0;
            [center, radius, metric]    = imfindcircles(me1,[12 30],'ObjectPolarity','dark'); % imfindcircles(abs(me1(1+B:end-B,1+B:end-B)./me2(1+B:end-B,1+B:end-B)),[12 30],'ObjectPolarity','bright');
        end
%             figure,imagesc(me1(1+B:end-B,1+B:end-B));axis image;colormap gray;viscircles(center(1,:), radius(1),'EdgeColor','g');
        center                      = round(center)+B;
        ROIr                        = center(1,2)-imSize/2:center(1,2)+imSize/2-1;
        ROIc                        = center(1,1)-imSize/2:center(1,1)+imSize/2-1;
        mov2                        = mov1(ROIr,ROIc,:);
        movOriginal                 = mov0(ROIr,ROIc,:);
    else
        ROIr                        = 510-imSize/2:510+imSize/2-1;
        ROIc                        = 401-imSize/2:401+imSize/2-1;
%         ROIr                        = 520-imSize/2:520+imSize/2-1;
%         ROIc                        = 410-imSize/2:410+imSize/2-1;
        mov2                        = mov1(ROIr,ROIc,:);
        movOriginal                 = mov0(ROIr,ROIc,:);
    end
  
%     me1                         = me;
%     me1(me1>0.40)              = 0.40;
%     for k = 1:settings.nFrames
% %         me2(:,:,k) = mov0(:,:,k)./me1;
%         subplot(221),imagesc(mov0(:,:,k)), axis image
%         subplot(222),imagesc(mov1(:,:,k)), axis image
%         subplot(223),imagesc(mov2(:,:,k)), axis image
% %         subplot(224),imagesc(me4(:,:,k)), axis image
%         drawnow,pause(0.03)
%     end

    %% Removing Static noise
    if 1
        [mov3,E,PP]                 = function_noiseRemoval(mov2, 20,1:20); % old 4um: (mov2, 20,[1:7,10:11])
        mov3                        = normalizingArray(mov3)-0.5;
    else
        mov3                        = normalizingArray(mov2);
    end
%     [movBG]                     = function_noiseRemoval_v2(movBG0, 3);
%     [movbgr,Er,PPr]           	= function_noiseRemoval_v2(movbgr0, 2);
%     [movbgc,Ec,PPc]           	= function_noiseRemoval_v2(movbgc0, 2);
    
%     figure
%     for k = 1:settings.nFrames
%         subplot(121),imagesc(mov2(:,:,k)),axis image
%         subplot(122),imagesc(mov3(:,:,k)),axis image,colormap(gray)
%         drawnow
%     end
    close all
    h = figure('Position',scrnsz);
    for k = 1:20
        subplot(5,4,k),imagesc(E(:,:,k)),axis image
    end
%     saveas(h,[FOLDER,'\',FILES{movnr,:}],'png')
    
    PPall(:,movnr) = PP;
end
save('PPall_4um','PPall','-v7.3')%,'p')


if 0
for movnr = 12
    %% FT filtering
    if 0
        L                           = 0.07;  % old 4um: L = 0.05     old 2um-3um: L = 0.15     2um-3um-4um: L = 0.1;                  11um: L = -0.05; 5um: L = -0.05 and -0.03;
        Lu                          = 0.915;  % 2um-3um-4um: Lu = 0.9;
        id_T                        = filterCustom(0,1,6,imSize,imSize); % old 2um-3um: filterCustom(0,1,6,60,60)	old 4um: filterCustom(0,2,7,60,60)  new 2um: filterCustom(0,1,8,60,60)	new 3um-4um: filterCustom(0,2,7,60,60)
%         figure
        for k = 1:settings.nFrames
            FRAME3                      = ifftshift(id_T).*fft2(mov3(:,:,k));
            frame30                     = ifft2(FRAME3);
            frame3                      = frame30;
%             frame3                      = normalizingArray(frame3);
            mov4(:,:,k)                 = frame3(:,:);
            frame3(frame3<L)            = L;
            frame3(frame3>Lu)           = Lu;
            mov5(:,:,k)                 = frame3(:,:);
%             subplot(131),imagesc(mov3(:,:,k)),axis image
%             subplot(132),imagesc(frame30),axis image%,colormap(gray)
%             subplot(133),imagesc(frame3),axis image%,colormap(gray)
%             drawnow,pause(0.03)
        end
    else
        mov4                        = mov3;
        mov5                        = mov3;
    end
    mov                         = normalizingArray(mov5);
    
    if 0
%         uncompressedVideo = VideoWriter('dropletTrappingAir.avi', 'Uncompressed AVI');
%         open(uncompressedVideo);
        fig = figure('position',[100 100 1200 600],'Color',[1.0 1.0 1.0]);
        for k = 30:size(mov,3)
            %             qualityCheck(:,:,movnr) = mov(:,:,k);
            %             subplot(5,6,movnr),imagesc(mov(:,:,k)),axis image
            %             subplot(121),
            A0                          = mov3(:,:,k);
            I(k)                        = mean(mean(mov3(700:800,700:800,1)))/mean(mean(A0(700:800,700:800)));
            A                           = A0;%*I(k);
%             A(800,800)                 	= 1;
            imagesc(A),axis image,colormap('gray')
            %             subplot(122),imagesc(mov2(:,:,k)),axis image,colormap('gray')
            drawnow
%             F = getframe(fig);
%             writeVideo(uncompressedVideo, F);
        end
%         close(uncompressedVideo);
    end
    
    return
    %% Viewing effect of code
%     figure
%     for k = 1%:settings.nFrames
%         figure(1), imagesc(mov0(ROIr,ROIc,k)), axis image%, title('Original frame')
%             set(gca,'xtick',[],'ytick',[]),colormap('gray')
%         figure(2), imagesc(mov3(ROIr,ROIc,k)), axis image%, title('Step 1: PCA filter')
%             set(gca,'xtick',[],'ytick',[]),colormap('gray')
%         figure(3), imagesc(mov4(:,:,k)), axis image%, title('Step 2: FT bandpass filter')
%             set(gca,'xtick',[],'ytick',[]),colormap('gray')
%         figure(4), imagesc(mov5(:,:,k)), axis image%, title('Step 3: Threshold filter')
%             set(gca,'xtick',[],'ytick',[]),colormap('gray')
%         figure(5), imagesc(E(ROIr,ROIc,1)), axis image%, title('Eigenvector 1')
%             set(gca,'xtick',[],'ytick',[]),colormap('gray')
%         figure(6), imagesc(E(ROIr,ROIc,2)), axis image%, title('Eigenvector 2')
%             set(gca,'xtick',[],'ytick',[]),colormap('gray')
%         figure(7), imagesc(E(ROIr,ROIc,3)), axis image%, title('Eigenvector 3')
%             set(gca,'xtick',[],'ytick',[]),colormap('gray')
%         figure(8), imagesc(E(ROIr,ROIc,4)), axis image%, title('Eigenvector 4+')
%             set(gca,'xtick',[],'ytick',[]),colormap('gray')
%         drawnow, pause(0.03)
%     end
% imagesc(E(:,:,5)), axis image, set(gca,'xtick',[],'ytick',[])
    %% Saving movie
    save(outputFile,'mov','movOriginal','settings','-v7.3')%,'p')
    close all
end
end