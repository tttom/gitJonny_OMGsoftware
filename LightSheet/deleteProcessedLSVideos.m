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
        folderNames='';
    end
    
    if (ischar(folderNames))
        folderNames={folderNames};
    end
    
    for (folderName=folderNames(:).')
        folderName=folderName{1};
        matFileList=dir(strcat(folderName,'/*.mat'));
        for (fileName={videoFileNameList.name})
            fileName=fileName{1}(1:end-4);
            filePathAndName=strcat(folderName,'/',fileName);
            configFileName=strcat(filePathAndName,'.json');
            videoFileName=strcat(filePathAndName,'.avi');
            disp(sprintf('Checking %s...',outputFileName))
            outputFileName=strcat(filePathAndName,'.mat');
            if (exist(videoFileName,'file') && exist(configFileName,'file'))
                storedVariables = whos('-file',outputFileName);
                if (ismember('recordedImageStack', {storedVariables.name}))
                    disp(sprintf('%s already processed.',outputFileName))
                    disp(sprintf('Deleting %s...',videoFileName))
                    delete(videoFileName);
                    disp(sprintf('%s Deleted.',videoFileName))
                end
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