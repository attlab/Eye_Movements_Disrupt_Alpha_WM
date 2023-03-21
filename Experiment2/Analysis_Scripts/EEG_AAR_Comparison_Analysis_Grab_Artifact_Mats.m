%{
EEG_AAR_Comparison_Analysis_Grab_Artifact_Mats
Author: Tom Bullock
Date: 10.14.22

Grab the artifact mats from the main EEG files to use in subsequent "no
AAR" analysis (while ensuring same trials get rejected).

%}

clear
close all

sourceDir = '/home/waldrop/Desktop/SCAMS/EEG_Prepro2_Avg_Baseline';
destDir = '/home/waldrop/Desktop/SCAMS/Data_Compiled';

subjects = [1:2,4:17,20,21:24,27:28]; % kicked out subjects 03 (bad EEG) and 19 (bad BEH)

for iSub=1:length(subjects)
    sjNum = subjects(iSub)
    for iSession =1:2
        load([sourceDir '/' sprintf('sj%02d_se%02d_wtfEye.mat',sjNum, iSession) ])
        all_arf(sjNum,iSession).arf = eeg.arf.artIndCleaned;
        all_arf(sjNum,iSession).bad_chans = eeg.badChannels;
        clear eeg
    end
end

save([destDir '/' 'EEG_Artifact_Mat.mat'],'all_arf','subjects')