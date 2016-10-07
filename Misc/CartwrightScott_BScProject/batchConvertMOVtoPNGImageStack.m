function batchConvertMOVtoPNGImageStack(inputFolder)

    %look in folder for all .mov files
    fileList=dir(strcat(inputFolder,'\*.mov'));
    
    %process each .mov file
    for fileIdx=1:length(fileList)
        fullFileName=strcat(inputFolder,'\',fileList(fileIdx).name);
        fprintf('Found .mov file: %s.\n',fullFileName);
        vidObj=VideoReader(fullFileName);
        
        %allocate memory for video frame
        nFrames=vidObj.NumberOfFrames;
        vidHeight=vidObj.Height;
        vidWidth=vidObj.Width;
        frame=zeros([vidHeight,vidWidth,3],'uint8');
        
        %create folder with same name as file
        mkdir(fullFileName(1:end-4));
        
        %load video frames into memory
        fprintf('Converting frames from video: %s...\n',fullFileName);
        for frameIdx=1:nFrames
            frame=read(vidObj,frameIdx);
            imwrite(frame,strcat(fullFileName(1:end-4),'/image_',num2str(100000+frameIdx),'.png'));
            fprintf('Converting frame %u of %u.\n',frameIdx,nFrames);
        end
        disp('Frames written.');
        
        %clear memory
        clear fullFileName vidObj nFrames vidHeight vidWidth frames
    end
end