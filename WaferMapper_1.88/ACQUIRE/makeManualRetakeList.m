%% Make manual retake list

%%Enter list
manualRetakeList = ...
[1:10:221];
%manualRetakeList = sort(manualRetakeList,'ascend')
TPN = GetMyDir;


save([TPN 'manualRetakeList.mat'],'manualRetakeList');