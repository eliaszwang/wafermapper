function [Files,MatFiles,labels]=GetSortedImagesAndMatfiles(directory)

filestruct = dir([directory filesep '*.tif']);
%labels=zeros(1,length(filestruct));
%Files=cell(1,length(filestruct));
%MatFiles=cell(1,length(filestruct));
  soNum = 0;
for i = 1:length(filestruct)
            %Extract Label
            
            Label = filestruct(i).name(length('SectionOverview_')+1:end-4);
          
            if  str2num(Label)>0
                soNum = soNum+1;
                Files{soNum}=[directory filesep filestruct(i).name];
                MatFiles{soNum}=[directory filesep filestruct(i).name(1:end-3) 'mat'];
                labels(soNum) = str2num(Label);
            end
end
[labels,indices]=sort(labels);
Files=Files(indices);
MatFiles=MatFiles(indices);