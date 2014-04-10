function[] = onePush(SPN,TPN)


clear pull mem


if ~exist(TPN,'dir') | ~exist(SPN,'dir')
    disp(error('missing directory'))
end

if ~exist(TPN,'dir')
    mkdir(targetFold)
end

memFileName = [SPN 'recordDataTransfer.mat'];

%% invarient Variables
delayCopySeconds = 30;  %Seconds between reading directories and starting file copy
delayCopyDays = delayCopySeconds/60/60/24;

%secondsBetweenFullCheck = 120;
%daysBetweenFullCheck = secondsBetweenFullCheck/60/60/24;




%% Read source folder

startFoldName = length(SPN) + 1;
startCheckTime =  datenum(clock);%record time at first read

%Update after first run through
currentCheckTime =  datenum(clock);
allGood = 1;
noNew = 1;

%% Read source folder
dirAnalysisStart = datenum(clock);
disp('Finding folders')
APN = findFolders(SPN); % find all folders in source directory

disp('Gathering File information')

pull.sourceDir = SPN;


pullSize = 0;
for f = 1: length(APN)  %run through every found folder
    sourceFolder = APN{f};
    targetFolder = [TPN sourceFolder(startFoldName:end)];
    pull.t(f).Folder = targetFolder;
    pull.f(f).Folder = sourceFolder;
    dAPN = dir(APN{f}); dAPN = dAPN(3:end);
    noDir = ~[dAPN.isdir];
    aFiles = dAPN(noDir);
    
    pull.f(f).Names = {aFiles.name}';
    pull.f(f).Sizes = [aFiles.bytes]';
    pullSize = pullSize + sum([aFiles.bytes]);
    pull.f(f).Times = [aFiles.datenum]';
    pull.f(f).lastChange = max([aFiles.datenum]);
    
    pull.t(f).Success = zeros(1,length(aFiles));
    pull.t(f).Sizes = zeros(1,length(aFiles));
    pull.t(f).Times = zeros(1,length(aFiles));
end

dirAnalysisStop = datenum(clock);
dirAnalysisTime = (dirAnalysisStop - dirAnalysisStart) * 24 * 60 * 60;
disp(['Analysis of ' SPN ' took ' num2str(dirAnalysisTime) ' seconds.'])

%% Get memory
if exist(memFileName,'file')
    load(memFileName)
    checkMemory = 1;
else
    checkMemory = 0;
    mem = pull;
end

%checkMemory =0


numAlreadyCopied = 0;
%%  Start Checking

totalMoved = 0;
rc = 0;  %Record number of popy reports
for f = 1:length(pull.f)
    %%Time out search
    
    if (length(pull.f(f).Times)) == (length(pull.f(f).Names))
        
        if ~exist(pull.t(f).Folder,'dir')
            mkdir(pull.t(f).Folder)
        end
        
        if exist('mem','var')  %if there is a remembered transfer
            memF = find(strcmp({mem.f.Folder}',pull.f(f).Folder)); %get folder id for mem
            
        else
            memF = []
        end
        
        %% Determine if anything in folder has been modified since last check
        folderChange = 1;
        if ~isempty(memF) & checkMemory
            if mem.f(memF).lastChange == pull.f(f).lastChange
                folderChange = 0;
            end
        end
        
        
        
        %% Does folder need to be checked
        if folderChange
            
            for i = 1 : length(pull.f(f).Names)
                
                
               
                
%                 if (datenum(clock) - pull.f(f).Times(i)) < delayCopyDays  %should wait for copy finish
%                     disp(sprintf('Copy of %s delayed until next round to prevent simultaneous read/write',dest))
%                     pause(delayCopyDays*24 * 60 * 60) 
%                 end
                source = [pull.f(f).Folder '\' pull.f(f).Names{i}];
                dest = [TPN source(startFoldName:end)];
                
                %% Check memory to see if file of right size and time has already been copied
                alreadyCopied = 0;
                 
                if checkMemory %%find remembered directory
                if ~isempty(memF) % if folder is remembered
                    if strcmp(mem.t(memF).Folder,pull.t(f).Folder); %if same target folder
                        memT = find(strcmp(mem.f(memF).Names,pull.f(f).Names(i)));  %find the name
                        if ~isempty(memT) %if file is remembered
                            %% Start selection
                            if mem.t(memF).Success(memT) %if previous copy was reported succesful
                                if mem.t(memF).Sizes(memT) == pull.f(f).Sizes(i)
                                    if mem.t(memF).Times(memT) >= pull.f(f).Times(i);
                                        alreadyCopied = 1;
                                    end
                                end
                            end
                        end
                    end
                end
                end
                
                if ~alreadyCopied
                    
                    
                    
                    %%Check destination to see if shouldCopy
                    shouldCopy = 0;
                    
                    if ~exist(dest,'file')
                        shouldCopy = 1;
                        disp(sprintf('create new file %s', dest))
                    else  %if file exists
                        fileInfo = dir(dest);
                        fileInfo = fileInfo(end);
                        
                        %%Check file info
                        if fileInfo.datenum< pull.f(f).Times(i) %replace old file
                            shouldCopy = 1;
                            disp(sprintf('replace old file %s', dest))
                        elseif (fileInfo.bytes ~= pull.f(f).Sizes(i)); %replace partial file
                            disp(sprintf('replaced partial file %s', dest))
                        else % if overwrite
                            %disp(sprintf('%s exists', dest))
                        end
                    end %whether dest file exists
                    
                    
                    if shouldCopy %copy if you should
                        status = 0;
                        noNew = 0;
                        startCopy = datenum(clock);
                        for c = 1:3
                            status = copyfile(source,dest); %copy the file
                            if status ~= 0  %make sure copy succeded
                                break
                            end
                            pause(1)
                        end %i status is bad
                        stopCopy = datenum(clock);
                        copySeconds = (stopCopy-startCopy)*24*60*60;
                        
                        if status ==0
                            disp('Failed to copy')
                            allGood = 0;
                            
                        else %if copy was a success
                            pull.t(f).Success(i) = 1;%if previous copy was reported succesful
                            pull.t(f).Sizes(i) = pull.f(f).Sizes(i);
                            pull.t(f).Times(i) = pull.f(f).Times(i);
                            if ( pull.f(f).Sizes(i)/1000000)<1
                                disp('Copy succeded')
                            else
                                disp(sprintf('%.1f MB copied at %.1f MB/sec',...
                                    pull.f(f).Sizes(i)/1000000, pull.f(f).Sizes(i)/1000000/copySeconds))
                            end
                        end
                        
                        
                    else
                        pull.t(f).Success(i) = 1;%if previous copy was reported succesful
                        pull.t(f).Sizes(i) = fileInfo.bytes;
                        pull.t(f).Times(i) = fileInfo.datenum;
                    end
                    
                else %if already copied
                    pull.t(f).Success(i) = 1;%if previous copy was reported succesful
                    pull.t(f).Sizes(i) = mem.t(memF).Sizes(memT);
                    pull.t(f).Times(i) =  mem.t(memF).Times(memT);
                    
                    numAlreadyCopied = numAlreadyCopied +1;
                    if ~mod(numAlreadyCopied,1000)
                        disp(sprintf(' %d files were already copied.',numAlreadyCopied))
                    end
                end
                
            totalMoved = totalMoved + pull.f(f).Sizes(i);  
            end %run all files in folder
            
                   %% save mem after each folder?
        if isempty(memF)
            mem.f(length(mem.f)+1) = pull.f(f);
            mem.t(length(mem.f)) = pull.t(f);
        else
            mem.f(memF) = pull.f(f); % Remember step.
            mem.t(memF) = pull.t(f);
        end
        safesave(memFileName,'mem');
            
            
        else %if all information is same length
                  totalMoved = totalMoved + sum(pull.f(f).Sizes);
        end

        
        
    end % same length
           disp(sprintf('%2.1f of %5.2f GB moved',(totalMoved/pullSize)*100,pullSize/10^9));    
 
end %run all folders

disp('Full check finished')
disp(sprintf(' %d files were already copied.',numAlreadyCopied))

if noNew
    disp(sprintf('No new files found after full check at %s', datestr(clock)))
else
    disp(sprintf('New files found after full check at %s', datestr(clock)))
end

%% finalize memory if function completes
mem = pull;
%safesave(memFileName,'mem');


