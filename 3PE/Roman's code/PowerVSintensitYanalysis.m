%Load the data
%clear all;
%workdir = 'Z:\Data\20150130\20xAir180nmSample\';
workdir = 'D:\Roman\Data\TF\ThreePhoton\20160505\NAD+\PowerScan3\';
files = sprintf('%s%s',workdir,'*tif');
%filelist =  dir('D:\Data\20150306\40xAir_WF_TF_9\*jpeg');
filelist =  dir(files);
s1 = size(filelist);
n_images = s1(1);
binning = 4;
%cropx = [400 700] / binning;
%cropy = [300 700] / binning;
cropx = [380 480] / binning;
cropy = [420 520] / binning;


clear profile;
clear slice;
clear imarray;
for i=1:n_images
    tf1 = imread(sprintf('%s%s',workdir,filelist(i).name));
    var = tf1; %% double(tf1(:,:,3))+double(tf1(:,:,2))+double(tf1(:,:,3));   
    imarray(:,:,i) = double(var(cropx(1):cropx(2),cropy(1):cropy(2)));
    profile(i) = sum(sum(var(cropx(1):cropx(2),cropy(1):cropy(2))))*0.001;
    %slice(i) = var(200,300);
    x_var(i) = i;
    %imarray(:,:,i) = double(var);% - double(bkg5sec2(:,:,1));
%     if i==round(n_images/2)
%         figure; imagesc(var);        
%     end
end
% for i=2:n_images
%     imarray(:,:,i) = imarray(:,:,i) - imarray(:,:,1);
% end
%norm_slice = normalise(slice-slice(n_images));
%%
cropx = [200 270];
cropy = [200 270];
clear imarray;
clear profile;
npoints = 13;
for i=1:npoints
    var = imread('F:\LightSheet\2016-05-11\NADH_2.3_790_1\NADH_2.3_790_1_MMStack.ome.tif',i);
    imarray(:,:,i) = double(var(cropx(1):cropx(2),cropy(1):cropy(2)));
    profile(i) = sum(sum(var(cropx(1):cropx(2),cropy(1):cropy(2))))*0.001;
end
profile = profile - profile(1);
%
figure; imagesc(var);
figure; imagesc(imarray(:,:,5));
figure; imagesc(imarray(:,:,10));
%profile = profile - profile(n_images); %% only if last image is a background !!!

%% Fit the data
%
y_rise = profile(2:npoints);
%x_power = 1200:-100:100;
x_power_crop =600:-100:100; 
y_rise = profile(npoints-5:npoints);
figure; plot(log(x_power),log(y_rise))
%y_fall = profile(2:npoints);
%hold on; plot(log(fliplr(x_power)),log(y_fall))
%%
%y_rise = profile(10:15);
%x_power = 800:100:1300;
%y_fall = profile(16:29);
figure; plot(log(x_power_crop),log(y_rise))
%hold on; plot(log(fliplr(x_power)),log(y_fall))
result = polyfit(log(x_power_crop),log(y_rise),1);
hold on; plot(log(x_power_crop),polyval(result,log(x_power_crop)));

%% Averaging over five measurements
n_average = 1;
n_points = n_images / n_average;
clear y_val;
clear y_std;
clear new_imarray;
for i=1:n_points
    y_val(i) = mean(profile((i-1)*n_average+1:i*n_average));
    y_std(i) = std(profile((i-1)*n_average+1:i*n_average))/sqrt(n_average);
    new_imarray(:,:,i) = mean(imarray(:,:,(i-1)*n_average+1:i*n_average),3);
end
%% Plot images in on figure
fig = figure;
n_points=13;
for i=1:n_points
    fig = subplot(3,5,i);
    imagesc(imarray(:,:,i));
    set(gca,'visible','off')
    axis equal;
end
colormap(gray);
