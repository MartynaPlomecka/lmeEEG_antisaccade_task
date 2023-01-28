srcDir = uigetdir;
files = dir(fullfile(srcDir,'*.mat'));
tab = [];
data = [];

for k = 1:length(files)
    baseFileName = files(k).name;
    fullFileName = fullfile(srcDir, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    load(fullFileName);
    
    %only stimlock, remove rows with saccades
    tmpinv=strcmp({stimEEG.event.type},'21') ;
    stimEEG.event(tmpinv)=[];
    clear tmpinv
    
    tmpinv=strcmp({stimEEG.event.type},'22') ;
    stimEEG.event(tmpinv)=[];
    clear tmpinv
    
    tmpinv=strcmp({stimEEG.event.type},'23') ;
    stimEEG.event(tmpinv)=[];
    clear tmpinv
    
    tmpinv=strcmp({stimEEG.event.type},'24') ;
    stimEEG.event(tmpinv)=[];
    clear tmpinv
    
    tmpinv=strcmp({stimEEG.event.cond},'') ;
    stimEEG.event(tmpinv)=[];
    stimEEG.data(:,:,tmpinv)=[];
    clear tmpinv
    
    tab = [tab;create_table(stimEEG)];
    data = cat(3,data,stimEEG.data);
    
end



function tab = create_table(EEG)

Condition = {EEG.event.cond}';
Age = ([EEG.event.age] .* 1)';
Dir = {EEG.event.dir}';

pieces = strsplit(EEG.fileName, '_');

ID = repelem(pieces(1), size(Age, 1))';

tab = table(Condition, ID, Age, Dir);

end