%{
EEG_Compute_Spectra
Author: Tom Bullock, UCSB Attention Lab
Date: 11.18.21

Compute spectra using FFTs

%}

function EEG_Compute_Spectra(sjNum)

% set dirs
%sourceDir = '/home/waldrop/Desktop/WTF_EYE/EEG_Prepro2_Avg_Baseline';
sourceDir = '/home/waldrop/Desktop/WTF_EYE/EEG_Prepro2_Avg_Baseline_noArtCorr';
%destDir = '/home/waldrop/Desktop/WTF_EYE/Spectra_Post_AAR';
destDir = '/home/waldrop/Desktop/WTF_EYE/Spectra_No_AAR';

% select subjects [MAKE THIS A FUNCTION WITH JOB SCRIPT] 
%sjNum = 1;% [1:7,9:14,16:20,22:27,31];


% merge data across sessions and split by condition
for iSession=1:2
    
    load([sourceDir '/' sprintf('sj%02d_se%02d_wtfEye.mat',sjNum,iSession)])
    
    % reject artifact trials from eeg and beh
    eegs = eeg.data;
    eegs = eegs(:,:,~eeg.arf.artIndCleaned);
    allBehStruct = allBehStruct(~eeg.arf.artIndCleaned);
    
    bothSessions(iSession).eegs=eegs;
    bothSessions(iSession).allBehStruct=allBehStruct;
    
    clear eegs allBehStruct
    
end

eegs = cat(3,bothSessions(1).eegs, bothSessions(2).eegs);
beh = cat(2,bothSessions(1).allBehStruct, bothSessions(2).allBehStruct);
clear bothSessions

% create condition vector
clear condBin
for i=1:length(beh)
    if      beh(i).memCondition==1 && beh(i).eyesCondition==2; condBin(i)=1;behBin(i)=beh(i);
    elseif  beh(i).memCondition==2 && beh(i).eyesCondition==2; condBin(i)=2;behBin(i)=beh(i);
    elseif  beh(i).memCondition==1 && beh(i).eyesCondition==6; condBin(i)=3;behBin(i)=beh(i);
    elseif  beh(i).memCondition==2 && beh(i).eyesCondition==6; condBin(i)=4;behBin(i)=beh(i);
    end
end

% split data by condition
for iCond=1:4
    allConds(iCond).eegs = eegs(:,:,find(condBin==iCond));
    allConds(iCond).beh = beh(find(condBin==iCond));
end


% get spectra using fft
for iCond=1:4
    
    % grab data
    eegData = allConds(iCond).eegs;
    
    % fft settings
    tStart = 1;
    tEnd = size(eegData,2);
    L = size(eegData,2);
    NFFT = L;
    sRate = 256;
    
    for iTrial = 1:size(eegData,3)
        for iChan = 1:size(eegData,1)
            spectra(iTrial,iChan,:) = fft(eegData(iChan,tStart:tEnd,iTrial),NFFT)/L;
        end
    end
    
    % get freqs
    freqs = 256/2*linspace(0,1,NFFT/2+1);
    
    % reduce freqs (do 4-20 Hz)
    spectra =spectra(:,:,find(freqs==1):find(freqs==30));
    
    freqs = freqs(find(freqs==1):find(freqs==30));
    
    % convert to power
    spectra = (abs(spectra)).^2;
    spectra = squeeze(mean(spectra,1));
    
    % compile averaged spectra into single matrix
    allSpectra(iCond,:,:) = spectra;
    
    clear eegData spectra 
    
end

save([destDir '/' sprintf('sj%02d_spectra.mat',sjNum)],'allSpectra','freqs')

return