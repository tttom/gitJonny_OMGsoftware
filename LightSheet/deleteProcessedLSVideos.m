%%% name:           deleteProcessedLSVideos
%%% author:         Jonathan Nylk
%%% date created:   14/08/2015
%%% description:    Deletes the large raw data .avi file associated with
%%%                 light-sheet images if the file has been processed and
%%%                 the data is already contained in the corresponding .mat
%%%                 file.
%%%                 The .avi file is replaced with a placeholder file.
%%%
%%% updates (latest first):
%%%
%%%
%%% END %%%

function deleteProcessedLSVideos(folderNames)

    if nargin<1
        folderNames={'F:\Stored Files\2015-08-11_Javier'};
    end
    
    if (ischar(folderNames))
        folderNames={folderNames};
    end
    
    for (folderName=folderNames(:).')
        folderName=folderName{1};
        disp(sprintf('Checking in folder: %s...',folderName))
        matFileList=dir(strcat(folderName,'/*.mat'));
        for (fileName={matFileList.name})
            fileName=fileName{1}(1:end-4);
            filePathAndName=strcat(folderName,'/',fileName);
            configFileName=strcat(filePathAndName,'.json');
            videoFileName=strcat(filePathAndName,'.avi');
            outputFileName=strcat(filePathAndName,'.mat');
            disp(sprintf('Checking %s...',outputFileName))
            if (exist(videoFileName,'file') && exist(configFileName,'file'))
                storedVariables = whos('-file',outputFileName);
                if (ismember('recordedImageStack', {storedVariables.name}))
                    [~,pos]=max(strcmp('recordedImageStack',{storedVariables.name}));
                    if (~isempty(storedVariables(5).size))
                        disp(sprintf('%s already processed.\nDeleting %s...',outputFileName))
%                         disp(sprintf('Deleting %s...',videoFileName))
                        delete(videoFileName);
                        placeHolderFileName=strcat(filePathAndName,'_PROCESSED.avi');
                        writerObj=VideoWriter(placeHolderFileName);
                        open(writerObj);
                        writeVideo(writerObj,zeros(2));
                        close(writerObj);
                        disp(sprintf('%s Deleted.',videoFileName))
                    else
                        disp(sprintf('"recordedImageStack" in %s is empty.\nNot deleting video file.',outputFileName))
                    end
                else
                    disp(sprintf('%s does not contain "recordedImageStack".\nNot deleting video file.',outputFileName))
                end
            else
                disp(sprintf('No matching video and/or config file for %s\nNot deleting video file.',outputFileName))
            end
        end
    end

    % Check for subfolders and handles these recursively
    directoryList=dir(folderName);
    for listIdx=1:length(directoryList)
        if directoryList(listIdx).isdir && directoryList(listIdx).name(1)~='.'
            expandedFolderName=strcat(folderName,'/',directoryList(listIdx).name);
            deleteProcessedLSVideos(expandedFolderName);
        end
    end
end