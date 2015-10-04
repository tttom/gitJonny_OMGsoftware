function quickStitchTwoMoviesMirrored(inputFolder1,inputFolder2,outputFolder)

    if nargin<1
        inputFolder1='F:\Stored Files\DavidLyons\Zebrafish2014-05-21\multicolor_7\Gauss_Rotation_HighRes_Raw';
    end
    if nargin<2
        inputFolder2='F:\Stored Files\DavidLyons\Zebrafish2014-05-21\multicolor_7\Airy_Rotation_HighRes_Raw';
    end
    if nargin<3
        outputFolder='F:\Stored Files\DavidLyons\Zebrafish2014-05-21\multicolor_7\AiryGaussMirrored_HighRes_Raw';
    end

    baseInputFileName1='Gauss_Rotation_HighRes';
    baseInputFileName2='Airy_Rotation_HighRes';
    baseOutputFileName='AiryGaussMirrored_HighRes';
    fileExt='.tif';
    for n=0:359 %360 deg rotation
        fileNo=num2str(10000+n);
        fileNum=fileNo(2:end);
        fullInputFileName1=strcat(inputFolder1,'\',baseInputFileName1,fileNum,fileExt);
        fullInputFileName2=strcat(inputFolder2,'\',baseInputFileName2,fileNum,fileExt);
        fullOutputFileName=strcat(outputFolder,'\',baseOutputFileName,fileNum,fileExt);
        
        inputImage1=imread(fullInputFileName1);
        inputImage2=fliplr(imread(fullInputFileName2));
        
        size1=size(inputImage1);
        size2=size(inputImage2);
        
        outputImage=zeros([size1(1),size1(2)+size2(2)],'uint8');
        outputImage(:,1:size1(2))=inputImage1;
        outputImage(:,size1(2)+1:end)=inputImage2;
        
%         imagesc(outputImage);axis image;colormap gray;drawnow;shg;
        imwrite(outputImage,fullOutputFileName);
        
    end
end