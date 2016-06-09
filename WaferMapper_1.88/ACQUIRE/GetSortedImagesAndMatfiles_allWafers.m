function [secOverviewFiles]=GetSortedImagesAndMatfiles_allWafers(waferList)

global GuiGlobalsStruct
UTSLdir = GuiGlobalsStruct.UTSLDirectory;

secOverviewFiles.labels = [];
secOverviewFiles.Files= {};
secOverviewFiles.MatFiles= {};
secOverviewFiles.sourceDir = {};
secOverviewFiles.alignedDir = {};

for w = 1:length(waferList)
wafDir = [UTSLdir filesep waferList{w} filesep];
soDir = [wafDir 'SectionOverviewsDirectory' filesep];
AlignSoDir = [wafDir  'SectionOverviewsAlignedWithTemplateDirectory' ];
if exist(AlignSoDir,'dir')
    newAlignSoDir = [AlignSoDir num2str(datenum(clock))];
    copyfile(AlignSoDir,newAlignSoDir)
end

filestruct = dir([soDir '*.tif']);
labels=zeros(1,length(filestruct));
Files=cell(1,length(filestruct));
MatFiles=cell(1,length(filestruct));
clear Files MatFiles labels sourceDir alignedDir

for i = 1:length(filestruct)
            %Extract Label
            Files{i}=[soDir filestruct(i).name];
            Label = filestruct(i).name(length('SectionOverview_')+1:end-4);
            MatFiles{i}=[soDir filestruct(i).name(1:end-3) 'mat'];
            labels(i) = str2num(Label);
            sourceDir{i} = soDir;
            alignedDir{i} = [AlignSoDir filesep];
            
end
[labels,indices]=sort(labels);
secOverviewFiles.labels = cat(2,secOverviewFiles.labels,labels);
secOverviewFiles.Files= cat(2,secOverviewFiles.Files,Files(indices));
secOverviewFiles.MatFiles= cat(2,secOverviewFiles.MatFiles,MatFiles(indices));
secOverviewFiles.sourceDir = cat(2,secOverviewFiles.sourceDir,sourceDir);
secOverviewFiles.alignedDir = cat(2,secOverviewFiles.alignedDir,alignedDir);
end