%{
Author: Tom Bullock
Date: 03.14.23

%}

clear
close all

sourceDir = '/home/waldrop/Desktop/SCAMS/BEH_Data_Processed';
destDir = '/home/waldrop/Desktop/SCAMS/Data_Compiled';

subjects = [1:2,4:17,20,21:24,27:28]; % kicked out subjects 03 (bad EEG) and 19 (bad BEH)

for iSub=1:length(subjects)
    
    sjNum = subjects(iSub)
    
    load([sourceDir '/' sprintf('sj%02d_allBeh.mat',sjNum)])
    
    for iCond=1:4
        for iSess=1:2
            if sjNum==20 && iCond==2 && iSess==1
                nBlocks=3;
            else
                nBlocks = 4;
            end
            n_bad_trials(iSess,iCond,iSub) =  length(masterStruct(iSess,iCond).allTrialData) -  (64*nBlocks);
        end
    end
    
end

% convert to proportion and get summary stats
m = squeeze(mean(n_bad_trials,1));
m = (m/256)*100;
props.m = m;
props.mean_bad_trials_prop = mean(m,2);
props.std_bad_trials_prop = std(m,0,2);
props.sem_bad_trials_prop = std(m,0,2)/sqrt(size(m,2));

% just get absolute numbers of trials (not proportions)
a = squeeze(sum(n_bad_trials,1));
counts.a = a;
counts.mean_bad_trials_n = mean(a,2);
counts.std_bad_trials_n = std(a,0,2);
counts.sem_bad_trials_n = std(a,0,2)/sqrt(size(a,2));



save([destDir '/' 'BAD_TRIALS_COUNT.mat'],'n_bad_trials','subjects','props','counts')




%mean_bad_trials_prop = mean_bad_trials/256*100;
%std_bad_trials_prop = std_bad_trials/