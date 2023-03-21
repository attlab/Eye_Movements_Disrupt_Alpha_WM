%{
EEG_Compute_Trial_Channel_Rej_Stats
Author:Tom Bullock
Date: 08.04.21

%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/EEG_Prepro2_Avg_Baseline';
destDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';

% set subjects
subjects = [1:7,9:14,16:20,22:27,31];

% loop through files
for iSub=1:length(subjects)
    sjNum = subjects(iSub)
    for iSession=1:2
        load([sourceDir '/' sprintf('sj%02d_se%02d_wtfEye.mat',sjNum,iSession)]);
        
        eyesCond = [allBehStruct.eyesCondition];
        memCond = [allBehStruct.memCondition];
        arfCond = eeg.arf.artIndCleaned;
        
        cnt1=0;cnt2=0;cnt3=0;cnt4=0;
        cnt5=0;cnt6=0;cnt7=0;cnt8=0;
        
        for iTrial=1:length(arfCond)
            
            if eyesCond(iTrial)==2 && memCond(iTrial)==1 && arfCond(iTrial)==1
                cnt1=cnt1+1;
            elseif eyesCond(iTrial)==6 && memCond(iTrial)==1 && arfCond(iTrial)==1
                cnt2=cnt2+1;
            elseif eyesCond(iTrial)==2 && memCond(iTrial)==2 && arfCond(iTrial)==1
                cnt3=cnt3+1;
            elseif eyesCond(iTrial)==6 && memCond(iTrial)==2 && arfCond(iTrial)==1
                cnt4=cnt4+1;
            end
            
            if eyesCond(iTrial)==2 && memCond(iTrial)==1
                cnt5=cnt5+1;
            elseif eyesCond(iTrial)==6 && memCond(iTrial)==1
                cnt6=cnt6+1;
            elseif eyesCond(iTrial)==2 && memCond(iTrial)==2
                cnt7=cnt7+1;
            elseif eyesCond(iTrial)==6 && memCond(iTrial)==2
                cnt8=cnt8+1;
            end
            
        end
        
        trialRejMat(iSub,iSession,:) = ([cnt1,cnt2,cnt3,cnt4]./[cnt5,cnt6,cnt7,cnt8])*100;
        
        % trialRejMat(iSub,iSession) =  (sum(eeg.arf.artIndCleaned)/length(eeg.arf.artIndCleaned))*100;
       
        % get channel rej stats
        chanRejMat(iSub,iSession) = length(eeg.badChannels);
        
    end
end

%meanTrialRejPercentageOverall = mean(mean(trialRejMat,1),2);

meanTrialRej = squeeze(mean(mean(trialRejMat,2),1));
semTrialRej = squeeze(std(mean(trialRejMat,2),0,1)./sqrt(25));

% get channel rej stats
meanChanRej = mean(mean(chanRejMat,1),2)/64*100;
semChanRej = std(mean(chanRejMat,2),0,1)./sqrt(25);

save([destDir '/' 'Rejection_Summary_Stats.mat'],'meanTrialRej','semTrialRej','meanChanRej','trialRejMat','chanRejMat')

