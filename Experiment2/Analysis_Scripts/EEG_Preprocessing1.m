%{
EEG_Preprocessing1
Author: Tom Bullock
Date: 11.12.19
Date last modified: 11.30.21
Import raw bdfs
Re-ref to mastoids
Add channel locations
Downsample to 256 Hz
Filter
Save as mat file
% REMEMBER TO BASELINE CORRECT IN SUBSEQUENT ANALYSIS SCRIPTS! %
%}

%function EEG_Preprocessing1(subjects)

clear
close all

% add paths
cd('/home/waldrop/Desktop/SCAMS/Analysis_Scripts');
cd('/home/waldrop/Desktop/SCAMS/Dependencies/eeglab14_1_2b');
eeglab
%clear
close all
cd('/home/waldrop/Desktop/SCAMS/Analysis_Scripts');

% choose subjects
subjects = [27]; 

% choose sessions
theseSessions = 2;%1:2; %[input 1:2 to do both, or 1 or 2 for separate]

% were data pre-merged? (subject,session on each row)
% preMergedDataSubjectsSessions = [11,1]; %[1,1;3,2];

% set dirs
rDir = '/home/waldrop/Desktop/SCAMS';
eegRawDir = [rDir '/' 'EEG_Raw'];
eegPrepro1Dir = [rDir '/' 'EEG_Prepro1'];
behDir = [rDir '/' 'BEH_Data_Processed'];

% loop through subs
for iSub=1:length(subjects)
    sjNum=subjects(iSub);
    
    % load behavioral data
    load([behDir '/' sprintf('sj%02d_allBeh.mat',sjNum)])
    
    % loop through sessions
    for iSession=theseSessions
        
%          % was data file pre-merged?
%          if ismember(sjNum,preMergedDataSubjectsSessions)
%             mergedFile=1;
%          else
%              mergedFile=0;
%          end
        
        % was data file pre-merged
        if sjNum==2 && iSession==1
            mergedFile=1;
        elseif sjNum==5 && iSession==2
            mergedFile=1;
        elseif sjNum==11 && iSession==1
            mergedFile=1;
        elseif sjNum==19 %%&& iSession==2
            mergedFile=1;
        elseif sjNum==20
            mergedFile=1;
        elseif sjNum==21
            mergedFile=1;
        elseif sjNum==27 && iSession==1
            mergedFile=1;
        else
            mergedFile=0;
        end
        
        % exception
        if sjNum==19 && iSession==1
            print('DO NOT RUN SJ19 SE01 AGAIN!!!!')
            return
            %break
        end
        
        if sjNum==11 && iSession==1
            print('DO NOT RUN SJ11 SE01 AGAIN!!!!')
            return
            %break
        end
        
        if sjNum==27 && iSession==1
            print('DO NOT RUN SJ27 SE01 AGAIN!!!!')
            return
            %break
        end
        
        % import\\load data
        if mergedFile
            load([eegRawDir '/' sprintf('sj%02d_se%02d_SCAMS.mat',sjNum,iSession) ]); 
        else
            EEG=pop_biosig([eegRawDir '/' sprintf('sj%02d_se%02d_SCAMS.bdf',sjNum,iSession)]); 
        end
        
        % downsample data
        EEG = pop_resample(EEG,256);
        
        % reference to mastoids
        EEG = pop_reref( EEG, [65 66],'keepref','on');
        
        % filter
        EEG = pop_eegfiltnew(EEG,0,80); % DO THIS FOR EOG
        
        % add channel locations
        EEG=pop_chanedit(EEG, 'lookup','/home/waldrop/Desktop/WTF_EYE/Dependancies/eeglab14_1_2b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp'); 
        
        % if merged file, convert EEG event codes from strings to nums
        
        if sjNum==19; mergedFile=0; end
        if sjNum==20 && iSession==1; mergedFile=0; end % exception becasue .mat but not merged
        if sjNum==21; mergedFile=0; end % exception because .mat but not merged
        
        if mergedFile
            for i=1:length(EEG.event)
                EEG.event(i).type = str2double(EEG.event(i).type);
            end
        end 
        
        % exception - for sj20 remove a block with no beh data
        if sjNum==20 && iSession==1
            EEG = pop_select(EEG,'nopoint',[EEG.event(4718).latency EEG.event(5122).latency-1]); 
            for i=1:length(EEG.event)
                EEG.event(i).type = str2double(EEG.event(i).type);
            end
        end
        
        if sjNum==20 && iSession==2
            EEG = pop_select(EEG,'nopoint',[EEG.event(5095).latency EEG.event(5443).latency-1]); 
            for i=1:length(EEG.event)
                EEG.event(i).type = str2double(EEG.event(i).type);
            end
        end
        
            
        % epoch files around fixation point (include ALL trials)
        EEG = pop_epoch(EEG,{102},[-.5 4]); % should output EEG!
        
        
        
        % prevent multiple targ triggers being coded into an epoch if trial
        % broken (which screws up later epoching around targets)
        for i=1:length(EEG.epoch)
            tmpVec = ismember(cell2mat(EEG.epoch(i).eventtype),[201,202,203,204,205,206,207,208]);
            %if tmpVec indicates more than one target code in epoch
            if sum(tmpVec)>1
                thisCode = EEG.epoch(i).eventtype(tmpVec);
                EEG.epoch(i).eventtype = thisCode(1); %broken trial to scrap all triggers except the first single target one
            end
        end
      
        % organize trial data in cronological order
        for iCond=1:4
            if iCond==1
                allBehStruct(1:size(masterStruct(iSession,cbOrder(iCond)).allTrialData,2)) = masterStruct(iSession,cbOrder(iCond)).allTrialData;
            else
                allBehStruct(length(allBehStruct)+1:length(allBehStruct)+length(masterStruct(iSession,cbOrder(iCond)).allTrialData))=masterStruct(iSession,cbOrder(iCond)).allTrialData;     
            end
        end
        
        % create broken trial vector
        cnt=0;
        for i=1:length(allBehStruct)          
           if allBehStruct(i).brokenTrial
               cnt=cnt+1;
              brokenTrialVec(cnt)=i;
           end    
        end
        
        %remove extra EEG trials if needed
%         EEG = pop_select(EEG,'notrial',215:228);
        
        % remove broken trials from EEG and trialData
        EEG = pop_select(EEG,'notrial',brokenTrialVec);
        allBehStruct(brokenTrialVec) = [];
        
        % re-epoch data to [-.5 to 2.5] around target onset
        EEG = pop_epoch(EEG,{201 202 203 204 205 206 207 208},[-.5 2.5]); % maybe extend out further pre-stim?
        
        
        % something strange going on with sj05se01 - broken trial vector
        % mismatch with EEG - for now, just kick out later trials
        if sjNum==5 && iSession==1
            EEG = pop_select(EEG,'notrial',728:897);
            allBehStruct(728:1024) = [];
        end
        
        % kick out one bad beh trial from sj28se01
        if sjNum==28 && iSession==1
            allBehStruct(737)=[];
        end
        
        
        % check that EEG and trial data match by checking event codes
        clear consistencyCheck
        for i=1:length(EEG.epoch)
           consistencyCheck(i,:) = [EEG.epoch(i).eventtype{1},allBehStruct(i).locTrigger,EEG.epoch(i).eventtype{1}-allBehStruct(i).locTrigger];          
        end
        

        
        if sum(consistencyCheck(:,3))~=0
            disp(sprintf('EEG_TRIAL DATA MISMATCH! INVESTIGATE sj%02d_se%02d!!',sjNum,iSession))
            return
        else
            disp('EEG\\TRIAL MATCH!! YAY!!')
        end
        
        % save the synced, epoched data
        save([eegPrepro1Dir '/' sprintf('sj%02d_se%02d_wtfEye_prepro1.mat',sjNum,iSession)],'EEG','allBehStruct','demographicInfo')
        
        clear brokenTrialVec  cnt consistencyCheck EEG thisCode tmpVec allBehStruct 
        
    end
end

%return