clc
clear

cd('\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Antisaccades\code\eeglab14_1_2b')
eeglab;
close all

x = dir('\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_new\data\EEG_preprocessed');
subjects = {x.name};
clear x

study_path = '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_new\';
analysis_path = [study_path 'scipts\data_processing_after_automagic\'];
cd(analysis_path)

subjects_exclude = [1,2];
subjects = subjects(~ismember(1:numel(subjects),subjects_exclude));
nsubjects = numel(subjects);

%%
for subj = 1:nsubjects
    
    
    datapath = strcat('\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_new\data\EEG_preprocessed\',subjects{subj});
    
    cd (datapath)
    
    if exist(strcat('gip_',subjects{subj},'_AS_EEG.mat')) > 0
        datafile= strcat('gip_',subjects{subj},'_AS_EEG.mat');
        load (datafile)
    elseif exist(strcat('oip_',subjects{subj},'_AS_EEG.mat')) > 0
        datafile= strcat('oip_',subjects{subj},'_AS_EEG.mat');
        load (datafile)
    else
        continue
    end
    
    
    %% REDUCE TO 105 electrodes
    
    tbl_channels = struct2table(EEG.chanlocs);
    el_excl = {'E48' 'E49' 'E56' 'E63' 'E68' 'E73' 'E81' 'E88' 'E94' 'E99' 'E107' 'E113' 'E119' 'E1' 'E8' 'E14' 'E17' 'E21' 'E25' 'E32' 'E125' 'E126' 'E127' 'E128'};
    
    
    tbl_excl = table(el_excl');
    tbl_excl.Properties.VariableNames = {'labels'};
    
    
    ind_excl = [];
    for i = 1:size(tbl_excl,1)
        ind_excl(end+1,1) = find(strcmp(tbl_excl.labels(i),[tbl_channels.labels]));
    end
    ind_excl';
    
    
    EEG = pop_select(EEG,'nochannel',ind_excl );
    
    %% Re-reference to average reference
    EEG = pop_reref(EEG,[]); %thats fine, its done AFTER AUTOMAGIC 
    
    %% unfold based cleaning
    mdl = [];
    mdl.toolboxpath =  '\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Antisaccades\unfold-develop-new';
    mdl.eventtypes = {'L_fixation'  'L_saccade'  'L_blink'};
    mdl.formula = {'y ~ 1'  'y ~ 1'  'y ~ 1'};
    mdl.categorical = {};
    mdl.codingschema = 'reference';
    mdl.timelimits = [-0.6000 1];
    mdl.amplitudeThreshold = 90;
    
    
    EEGorg = EEG;
    
    addpath(mdl.toolboxpath);
    init_unfold; %initialize the toolbox
    
    cfg = []; % Model specification
    
    cfg.eventtypes = mdl.eventtypes;
    
    cfg.formula = mdl.formula;
    
    cfg.categorical = mdl.categorical;
    
    cfg.codingschema = mdl.codingschema;
    
    EEG = uf_designmat(EEG,cfg); % create design matrix
    
    cfg = []; % Time expansion
    
    cfg.timelimits = mdl.timelimits;
    
    EEG = uf_timeexpandDesignmat(EEG,cfg); % time-expand the design matrix
    
    %% Fit the model
    
    EEG = uf_glmfit(EEG); % Model fitting, solve the linear model (with deconvolution)
    ufresult = uf_condense(EEG);
    
    %% Make cleaning-model-predicted and cleaned continuous EEG
    
    desMat = ufresult.unfold.Xdc; %[cntSp x (epoSp * param)], complete design matrix (2D)
    
    beta = permute(ufresult.beta,[2,3,1]); %[epoSp x params x chan], model parameters (need to reorder dimensions for matrix multiplication later on)
    
    % beta = beta(:,indKeepBeta,:); %]epoSp x params x chan], only keep model parameters of interest
    
    beta = reshape(beta,[],EEG.nbchan); %[(epoSp * param) x chan], make model parameters 2D for matrix multiplication
    
    cntModel = (desMat * beta)'; % [chan x cntSp], matrix multiplication and reorientation, this is the predicted timeseries
    
    EEG_model = EEGorg; %model-predicted EEG
    
    EEG_model.data = cntModel;
    
    EEG_res = EEGorg; %residual EEG, what is left when removing the model prediction
    
    EEG_res.data = EEGorg.data - EEG_model.data;
    
    EEG_res = EEG;
    
    
    %% triggers renaming
    countblocks = 1;
    for e = 1:length(EEG.event)
        if strcmp(EEG.event(e).type, 'boundary')
            countblocks = countblocks + 1;
            continue;
        end
        if countblocks == 2 || countblocks == 3 || countblocks == 4 % antisaccade blocks
            if strcmp(EEG.event(e).type,'10  ') % change 10 to 12 for AS
                EEG.event(e).type = '13';
            elseif strcmp(EEG.event(e).type,'11  ')
                EEG.event(e).type = '14'; % change 11 to 13 for AS
            end
            
            if strcmp(EEG.event(e).type,'40  ')
                EEG.event(e).type = '41  ';
            end
            
        end
        
        if countblocks == 1 || countblocks == 5 %prosaccade blocks
            if strcmp(EEG.event(e).type,'10  ') % change numbers to words
                EEG.event(e).type = '11';
            elseif strcmp(EEG.event(e).type,'11  ')
                EEG.event(e).type = '12'; % change numbers to words
            end
        end
        
        
    end
    
    EEG.event(strcmp('boundary',{EEG.event.type})) = [];
    rmEventsIx = strcmp('L_fixation',{EEG.event.type});
    rmEv =  EEG.event(rmEventsIx);
    EEG.event(rmEventsIx) = [];
    EEG.event(1).dir = []; %left or right
    EEG.event(1).cond = [];%pro or anti
    %% rename EEG.event.type
    previous = '';
    for e = 1:length(EEG.event)
        if strcmp(EEG.event(e).type, 'L_saccade')
            if strcmp(previous, '11')
                EEG.event(e).type = '21';
                EEG.event(e).cond = 'pro';
                EEG.event(e).dir = 'left';
                
                %pro left
            elseif strcmp(previous, '12')
                EEG.event(e).type = '22';
                EEG.event(e).cond = 'pro';
                EEG.event(e).dir = 'right';
                
            elseif strcmp(previous, '13')
                EEG.event(e).type = '23';
                EEG.event(e).cond = 'anti';
                EEG.event(e).dir = 'left';
                
            elseif strcmp(previous, '14')
                EEG.event(e).type = '24';
                EEG.event(e).cond = 'anti';
                EEG.event(e).dir = 'right';
                
            else
                EEG.event(e).type = 'invalid';
            end
            
        end
        if ~strcmp(EEG.event(e).type, 'L_fixation') ...
                && ~strcmp(EEG.event(e).type, 'L_blink')
            previous = EEG.event(e).type;
        end
    end
    
    %% remove everything from EEG.event which is not saccade or trigger sent by me from the stimulus pc
    
    tmpinv=strcmp({EEG.event.type}, 'invalid') | strcmp({EEG.event.type}, 'L_blink') ;
    EEG.event(tmpinv)=[];
    
    %% renaming errors
    
    
    for e = 1:length(EEG.event)
        
        if strcmp(EEG.event(e).type, '23') && ( EEG.event(e).sac_startpos_x > EEG.event(e).sac_endpos_x)
            EEG.event(e).accuracy = 'error_anti_sacc';
            EEG.event(e-1).accuracy = 'error_anti_cue';
            
        elseif strcmp(EEG.event(e).type, '23') && (EEG.event(e).sac_startpos_x < EEG.event(e).sac_endpos_x)
            EEG.event(e).accuracy = 'correct_anti_sacc';
            EEG.event(e-1).accuracy = 'correct_anti_cue';
            
        elseif strcmp(EEG.event(e).type, '24') && (EEG.event(e).sac_startpos_x <EEG.event(e).sac_endpos_x)
            EEG.event(e).accuracy = 'error_anti_sacc';
            EEG.event(e-1).accuracy = 'error_anti_cue';
            
        elseif strcmp(EEG.event(e).type, '24') && (EEG.event(e).sac_startpos_x >EEG.event(e).sac_endpos_x)
            EEG.event(e).accuracy = 'correct_anti_sacc';
            EEG.event(e-1).accuracy = 'correct_anti_cue';
            
        elseif strcmp(EEG.event(e).type, '21') && ( EEG.event(e).sac_startpos_x < EEG.event(e).sac_endpos_x)
            EEG.event(e).accuracy = 'error_pro_sacc';
            EEG.event(e-1).accuracy = 'error_pro_cue';
            
        elseif strcmp(EEG.event(e).type, '21') && (EEG.event(e).sac_startpos_x > EEG.event(e).sac_endpos_x)
            EEG.event(e).accuracy = 'correct_pro_sacc';
            EEG.event(e-1).accuracy = 'correct_pro_cue';
            
        elseif strcmp(EEG.event(e).type, '22') && (EEG.event(e).sac_startpos_x >EEG.event(e).sac_endpos_x)
            EEG.event(e).accuracy = 'error_pro_sacc';
            EEG.event(e-1).accuracy = 'error_pro_cue';
            
        elseif strcmp(EEG.event(e).type, '22') && (EEG.event(e).sac_startpos_x <EEG.event(e).sac_endpos_x)
            EEG.event(e).accuracy = 'correct_pro_sacc';
            EEG.event(e-1).accuracy = 'correct_pro_cue';
            
            
        else
            EEG.event(e).accuracy = 'NA';
        end
    end
    
    %% CALCULATE RT FOR EACH TRIALS
    
    
    for e = 1:length(EEG.event)
        
        if strcmp(EEG.event(e).type, '23')
            EEG.event(e).rt = (EEG.event(e).latency - EEG.event(e-1).latency)*2;%for sacc
            EEG.event(e-1).rt = (EEG.event(e).latency - EEG.event(e-1).latency)*2; %for the "pair"cue
            
            
        elseif strcmp(EEG.event(e).type, '24')
            EEG.event(e).rt = (EEG.event(e).latency - EEG.event(e-1).latency)*2;
            EEG.event(e-1).rt = (EEG.event(e).latency - EEG.event(e-1).latency)*2;
            
            
        elseif strcmp(EEG.event(e).type, '21')
            EEG.event(e).rt = (EEG.event(e).latency - EEG.event(e-1).latency)*2;
            EEG.event(e-1).rt = (EEG.event(e).latency - EEG.event(e-1).latency)*2;
            
            
        elseif strcmp(EEG.event(e).type, '22')
            EEG.event(e).rt = (EEG.event(e).latency - EEG.event(e-1).latency)*2;
            EEG.event(e-1).rt = (EEG.event(e).latency - EEG.event(e-1).latency)*2;
            
            
        else
            EEG.event(e).rt = 'NA';
        end
    end
    
    %% amplitude too small
    tmperrsacc6=find(strcmp({EEG.event.type}, '22') ...
        & [EEG.event.sac_amplitude]<1.5);
    tmperrsacc7=find(strcmp({EEG.event.type}, '21') ...
        & [EEG.event.sac_amplitude]<1.5);
    tmperrsacc8=find(strcmp({EEG.event.type}, '23') ...
        & [EEG.event.sac_amplitude]<1.5);
    tmperrsacc9=find(strcmp({EEG.event.type}, '24') ...
        & [EEG.event.sac_amplitude]<1.5);
    tmperr69=[tmperrsacc6 (tmperrsacc6-1) tmperrsacc7 (tmperrsacc7-1) tmperrsacc8 (tmperrsacc8-1) tmperrsacc9 (tmperrsacc9-1)];
    EEG.event(tmperr69)=[];
    
    clear tmperrsacc1 tmperrsacc2 tmperrsacc3 tmperrsacc4 tmperrsacc6 tmperrsacc7 tmperrsacc8 tmperrsacc9
    
    
    
    %% delete cues where there was no saccade afterwards
    
    tmperrcue10=  find(strcmp({EEG.event.type}, '11')) ;
    for i=1:length(tmperrcue10)
        pos = tmperrcue10(i);
        if ~ (strcmp(EEG.event(pos+1).type , '21'))
            
            EEG.event(pos).type='missingsacc'; %cue
        end
    end
    
    %%11
    tmperrcue11 =   find(strcmp({EEG.event.type}, '12'))    ;
    for i=1:length(tmperrcue11)
        pos = tmperrcue11(i);
        if ~ (strcmp(EEG.event(pos+1).type , '22'))
            
            EEG.event(pos).type='missingsacc'; %cue
        end
    end
    
    
    tmperrcue12=  find(strcmp({EEG.event.type}, '13')) ;
    for i=1:length(tmperrcue12)
        pos = tmperrcue12(i);
        if ~ (strcmp(EEG.event(pos+1).type , '23'))
            
            EEG.event(pos).type='missingsacc'; %cue
        end
    end
    
    %%11
    tmperrcue13 =   find(strcmp({EEG.event.type}, '14'))    ;
    for i=1:length(tmperrcue13)
        pos = tmperrcue13(i);
        if ~ (strcmp(EEG.event(pos+1).type , '24'))
            
            EEG.event(pos).type='missingsacc'; %cue
        end
    end
    
    tmpinv=strcmp({EEG.event.type}, 'missingsacc') ;
    EEG.event(tmpinv)=[];
    
    
    %% delete saccades and cues when the saccade comes faster than 100ms after cue
    tmpevent=length(EEG.event);
    saccpro=find(strcmp({EEG.event.type},'22')==1 | strcmp({EEG.event.type},'21')==1); % find rows where there is a saccade
    saccanti=find(strcmp({EEG.event.type},'24')==1 | strcmp({EEG.event.type},'23')==1);%find rows where there is a saccade
    
    for b=1:size(saccpro,2)
        
        if (EEG.event(saccpro(1,b)).latency-EEG.event(saccpro(1,b)-1).latency)<100 %50 because 100ms
            EEG.event(saccpro(b)).type='too_fast'; %saccade
            EEG.event(saccpro(b)-1).type = 'too_fast'; %cue
        end
    end
    
    for b=1:size(saccanti,2)
        
        if (EEG.event(saccanti(b)).latency-EEG.event(saccanti(1,b)-1).latency)<100
            EEG.event(saccanti(b)-1).type ='too_fast';
            EEG.event(saccanti(b)).type ='too_fast';
        end
        
    end
    
    tmpinv=strcmp({EEG.event.type}, 'too_fast') ;
    EEG.event(tmpinv)=[];
    clear tmpinv
    
    %% remove everything except stim-response
    tmpinv=strcmp({EEG.event.rt}, 'NA') ;
    EEG.event(tmpinv)=[];
    clear tmpinv
    
    %% add trial number
    
    n =length(EEG.event)/2;
    a = reshape( repmat( [1:n], 2,1 ), 1, [] );
    a = a';
    
    for i =1:length(EEG.event)
        EEG.event(i).trial = a(i);
    end
    
    %% add age info
    id = regexp(EEG.comments(1,:), '.*All_Subjects[\/\\](?<ID>.*)[\/\\].*', 'names').ID;
    is_old = any(ismember(OLD.Subject, id));
    for e=1:size(EEG.event,2)
        EEG.event(e).age = is_old;
    end
    
   
    %%
    tmpinv=strcmp({EEG.event.accuracy}, 'error_pro_cue');
    EEG.event(tmpinv)=[];
    
    tmpinv=strcmp({EEG.event.accuracy}, 'error_pro_sacc');
    EEG.event(tmpinv)=[];
    
    tmpinv=strcmp({EEG.event.accuracy}, 'error_anti_cue');
    EEG.event(tmpinv)=[];
    
    tmpinv=strcmp({EEG.event.accuracy}, 'error_anti_sacc');
    EEG.event(tmpinv)=[];
    %%
    saccEEG = pop_epoch( EEG, {'22','21','23','24'}, [-0.6  0.6], 'epochinfo', 'yes');
    saccEEG = pop_rmbase(saccEEG,[-200 0]);
    ERPsacc(:,:) = mean(saccEEG.data(:,:,:),3);
   
    stimEEG = pop_epoch( EEG, {'12','11','13','14'}, [-0.6  0.6], 'epochinfo', 'yes');
    stimEEG = pop_rmbase(stimEEG,[-200 0]);
    ERPstim(:,:) = mean(stimEEG.data(:,:,:),3);
    %% save this ready to statistics data structure
    
      mkdir('\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_new\data\mass_univ_tables\segmented_data_correcttrials_only_baselineremov\', id)
      save(['\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_new\data\mass_univ_tables\segmented_data_correcttrials_only_baselineremov\' id '\' id '_sacclockedEEG.mat'], 'saccEEG');
      save(['\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_new\data\mass_univ_tables\segmented_data_correcttrials_only_baselineremov\' id '\' id '_stimlockedEEG.mat'], 'stimEEG');
      save(['\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_new\data\mass_univ_tables\segmented_data_correcttrials_only_baselineremov\' id '\' id '_stimlockedERP.mat'], 'ERPstim');
      save(['\\psyger-stor02.d.uzh.ch\methlab\Neurometric\Anti_new\data\mass_univ_tables\segmented_data_correcttrials_only_baselineremov\' id '\' id '_sacclockedERP.mat'], 'ERPsacc');

end

