function SpatialIEM_TrainBoth_Balance_4conds (subjects)
%==========================================================================
%{
Purpose: run spatial encoding model on single bandpass e.g. 8-12 Hz

Original Author:
Joshua J. Foster
joshua.james.foster@gmail.com
University of Chicago
August 12, 2015
Modified by Tom Bullock
UCSB Attention Lab
%}
%==========================================================================

close all

% quick overwrite
%subjects = [1:10,13,16,19:24];


subs = subjects;


nSubs = length(subs);

%% process for IEM (1) or HILBERT (2)
processForIEM=1;

% setup directories
%root = '/home/garrett/WTF_Bike';
root = '/home/waldrop/Desktop/WTF_EYE';
%dRoot = [root '/Data/'];
dRoot = [root '/BEH'];

% % setup directories
%
% %root = pwd;
% out = 'AnalysisScripts';
% dRoot = [root(1:end-length(out)),'Data/'];

if processForIEM==1
    %eRoot = [root(1:end-length(out)),'EEG/'];
    %eRoot = [root '/EEG/'];
    eRoot = [root '/home/waldrop/Desktop/WTF_EYE/EEG_Prepro2_Avg_Baseline/'];
else
    eRoot = [root '/HILBERT/'];
end

%hRoot = [root(1:end-length(out)),'HILBERT/'];
%hRoot = [root '/HILBERT/'];

%%%CLEAN UP%%%

% choose method for finding the min number of trials per location bin
% (1=per-condition, 2=across conditions)
findMinType = 1; %********************************JORDAN EDIT, was 1 now 2**************

% % load matrix of min trials per loc bin per condition
% if processForIEM == 1
%     minLocsData = load([root '/Analysis_Scripts/Minimum_Location_Bin_Mat_AccTrials.mat']);
% elseif processForIEM== 2
%     minLocsData = load([root '/Analysis_Scripts/HILBERT_Minimum_Location_Bin_Mat_AccTrials.mat']);
% end

for bandpassLoop=1;%bands % loop through different bandpass analyses e.g. alpha, theta
    
    if bandpassLoop==1
        name = '_SpatialTF_ALPHA.mat'; % name of files to be saved for IEM
        hilbName = '_Hilbert_ALPHA.mat';
        thisFreq = {'Alpha'};
        freqBandpass = [8 12];
    else
        name = '_SpatialTF_THETA.mat';
        hilbName = '_Hilbert_THETA.mat';
        thisFreq = {'Theta'};
        freqBandpass = [4 7];
    end
    
    % Loop through participants
    %matlabpool open 72
    for s = 1:nSubs
        
        % had to move all the settings insdie the parfor to get it to work...
        em = [];
        EEG = [];
        eegInfo = [];
        WeightsTotal = [];
        allWeightsTotal = [];
        
        % parameters to set
        em.nChans = 8; % # of channels;
        em.nBins = em.nChans; % # of stimulus bins, was em.nChans; using 16 since doing both conditions. > 8 is exercise data
        em.nIter = 10; % # of iterations %%% WAS SET TO 10!!!
        em.nBlocks = 8;%4; %3; % # of blocks for cross-validation
        em.frequencies = freqBandpass; % frequency bands to analyze
        em.bands = thisFreq;
        em.window = 4;
        em.Fs = 250; % WAS 250
        em.nElectrodes = 999; % get later
        em.time = -.5*1000:1000/em.Fs:1.9961*1000; %   -500:4:2000; % time points of interest
        
        % for brevity in analysis
        nChans = em.nChans;
        nBins = em.nBins;
        nIter = em.nIter;
        nBlocks = em.nBlocks;
        freqs = em.frequencies;
        times = em.time;
        nFreqs = size(em.frequencies,1);
        nElectrodes = em.nElectrodes;
        nSamps = length(em.time);
        Fs = em.Fs;
        
        % Specify basis set
        em.sinPower = 7;
        em.x = linspace(0, 2*pi-2*pi/nBins, nBins);
        em.cCenters = linspace(0, 2*pi-2*pi/nChans, nChans);
        em.cCenters = rad2deg(em.cCenters);
        pred = sin(0.5*em.x).^em.sinPower; % hypothetical channel responses
        pred = wshift('1D',pred,5); % shift the initial basis function
        basisSet = nan(nChans,nBins);
        for c = 1:nChans
            basisSet(c,:) = wshift('1D',pred,-c); % generate circularly shifted basis functions
        end
        em.basisSet = basisSet; % save basis set to data structure
        
        sn = subs(s);
        fprintf('Subject:\t%d\n',sn)
        
        % Grab data------------------------------------------------------------
        
        
        % TOM ADDED FROM HERE ON
        eegDir = '/home/waldrop/Desktop/WTF_EYE/EEG_Prepro2_Avg_Baseline';
        
        % merge data across sessions and split by condition
        for iSession=1:2
            
            % load data
            load([eegDir '/' sprintf('sj%02d_se%02d_wtfEye.mat',sn,iSession)])
            
            % reject artifact trials from eeg and beh
            eegs = eeg.data;
            eegs = eegs(:,:,~eeg.arf.artIndCleaned);
            allBehStruct = allBehStruct(~eeg.arf.artIndCleaned);
            
            % combine sessions
            bothSessions(iSession).eegs=eegs;
            bothSessions(iSession).allBehStruct=allBehStruct;
            
            clear eegs allBehStruct
            
        end
        
        % create merged session matrices
        eegs = cat(3,bothSessions(1).eegs, bothSessions(2).eegs); % [elects x times x trials]
        beh = cat(2,bothSessions(1).allBehStruct, bothSessions(2).allBehStruct); % [trials]
        clear bothSessions
        
        % create vector of condition labels
        clear condBin
        for i=1:length(beh)
            if      beh(i).memCondition==1 && beh(i).eyesCondition==2; condBin(i)=1;%behBin(i)=beh(i);
            elseif  beh(i).memCondition==2 && beh(i).eyesCondition==2; condBin(i)=2;%behBin(i)=beh(i);
            elseif  beh(i).memCondition==1 && beh(i).eyesCondition==6; condBin(i)=3;%behBin(i)=beh(i);
            elseif  beh(i).memCondition==2 && beh(i).eyesCondition==6; condBin(i)=4;%behBin(i)=beh(i);
            end
        end
        
        % split eeg and beh data by condition
        for iCond=1:4
            allConds(iCond).eegs = eegs(:,:,find(condBin==iCond));
            allConds(iCond).beh = beh(find(condBin==iCond));
        end
        
        clear condBin
        
        % find location bin with min count per condition (to balance IEM)
        for iCond=1:4
            
            clear beh posBin minCnt
            beh=allConds(iCond).beh;
            
            for i=1:length(beh)
                posBin(i)=beh(i).stimLoc;
            end
            
            binCnt = [];
            for bin = 1:8
                binCnt(bin) = sum(posBin == bin);
            end
            
            allConds(iCond).minCnt = min(binCnt); % add minCnt to allConds structure
            allPosBin(iCond).posBin=posBin; % create structure with just position bins
            
            clear posBin
        end
        
        % find min bin cnt across conditions (-1 ensures variable
        minCnt = min([allConds.minCnt]) - 1;
        
        
        
%         % Get position bin index from behavior file
%         restfName = [dRoot, sprintf('sj%02d_exerCon01_changeDect_MixModel_wBias_accTrials.mat',subs(s))];
%         exfName = [dRoot, sprintf('sj%02d_exerCon02_changeDect_MixModel_wBias_accTrials.mat',subs(s))];

%         cond1fName = [dRoot, sprintf('/sj%d01_newBeh',sn)];
%         cond2fName = [dRoot, sprintf('/sj%d02_newBeh',sn)];
%         cond3fName = [dRoot, sprintf('/sj%d03_newBeh',sn)];
%         cond4fName = [dRoot, sprintf('/sj%d04_newBeh',sn)];
        
        % load files
%         cond1_tmp = []; cond1_beh = [];
%         cond1_tmp = load(cond1fName);
        %cond1_beh = [cond1_tmp.trialInfoUnbroken.stimLoc];
        cond1_beh = [allConds(1).beh.stimLoc];
        
%         cond2_tmp = []; cond2_beh = [];
%         cond2_tmp = load(cond2fName);
        %cond2_beh = [cond2_tmp.trialInfoUnbroken.stimLoc];
        cond2_beh = [allConds(2).beh.stimLoc];
        
%         cond3_tmp = []; cond3_beh = [];
%         cond3_tmp = load(cond3fName);
        %cond3_beh = [cond3_tmp.trialInfoUnbroken.stimLoc];
        cond3_beh = [allConds(3).beh.stimLoc];
        
%         cond4_tmp = []; cond4_beh = [];
%         cond4_tmp = load(cond4fName);
        %cond4_beh = [cond4_tmp.trialInfoUnbroken.stimLoc];
        cond4_beh = [allConds(4).beh.stimLoc];
        
        
        %beh.trial.posBin = [rest_beh.trial.posBin,ex_beh.trial.posBin]; %merge data
        clear beh
        beh.trial.posBin = [cond1_beh,cond2_beh,cond3_beh,cond4_beh];
        
        %Jordan Added: create index to keep track of rest and exercise position trials
        %trial_indx = [repmat(1,1,length(rest_beh.trial.posBin)), repmat(2,1,length(ex_beh.trial.posBin))];
        trial_indx = [repmat(1,1,length(cond1_beh)),repmat(2,1,length(cond2_beh)),repmat(3,1,length(cond3_beh)),repmat(4,1,length(cond4_beh))];
        
        em.posBin = beh.trial.posBin'; % add to fm structure so it's saved
        clear posBin
        posBin = em.posBin;
        
        % Get EEG data
        %restfName = [eRoot, sprintf('sj%02d_exerCon01_changeDect_EEG_accTrials.mat',subs(s))];
        %exfName = [eRoot, sprintf('sj%02d_exerCon02_changeDect_EEG_accTrials.mat',subs(s))];
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%5
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        %%%HEREE%%%
        
%         cond1fName = [eRoot, sprintf('%d01_EEG.mat',subs(s))];
%         cond2fName = [eRoot, sprintf('%d02_EEG.mat',subs(s))];
%         cond3fName = [eRoot, sprintf('%d03_EEG.mat',subs(s))];
%         cond4fName = [eRoot, sprintf('%d04_EEG.mat',subs(s))];
        
        % load file
%         cond1_tmp = []; cond1_eeg = [];
%         cond1_tmp = load(cond1fName);
        %cond1_eeg = cond1_tmp.eeg;
        cond1_eeg = allConds(1).eegs; % chans x times x trials
        
%         cond2_tmp = []; cond2_eeg = [];
%         cond2_tmp = load(cond2fName);
        %cond2_eeg = cond2_tmp.eeg;
        cond2_eeg = allConds(2).eegs; % chans x times x trials
        
%         cond3_tmp = []; cond3_eeg = [];
%         cond3_tmp = load(cond3fName);
        %cond3_eeg = cond3_tmp.eeg;
        cond3_eeg = allConds(3).eegs; % chans x times x trials
        
%         cond4_tmp = []; cond4_eeg = [];
%         cond4_tmp = load(cond4fName);
        %cond4_eeg = cond4_tmp.eeg;
        cond4_eeg = allConds(4).eegs; % chans x times x trials
        
        
        % permute to correct format (trials x chans x times)
        cond1_eeg = permute(cond1_eeg,[3,1,2]);
        cond2_eeg = permute(cond2_eeg,[3,1,2]);
        cond3_eeg = permute(cond3_eeg,[3,1,2]);
        cond4_eeg = permute(cond4_eeg,[3,1,2]);
        %%%%START%%%
        
        
        
        
        
        
        
        
        %%%HERE%%%
        
        
        %Jordan Added
        %eeg.data = cat(1,cond1_eeg.data,cond2_eeg.data,cond3_eeg.data,cond4_eeg.data); %merge eeg trials along first dimension (trials x channels x timepoints)
        
        eeg.preTime = eeg.preTime;
        eeg.postTime = eeg.postTime;
        eeg.chanLabels = eeg.chanLabels;
        eeg.sampRate = eeg.sampRate;
        
        eeg.data = [];
        eeg.data = cat(1,cond1_eeg,cond2_eeg,cond3_eeg,cond4_eeg);
        %eeg.arf.artIndCleaned = [cond1_eeg.arf.artIndCleaned, cond2_eeg.arf.artIndCleaned, cond3_eeg.arf.artIndCleaned, cond4_eeg.arf.artIndCleaned];
        
%         eeg.preTime = cond1_eeg.preTime;
%         eeg.postTime = cond1_eeg.postTime;
%         eeg.chanLabels = cond1_eeg.chanLabels;
%         eeg.sampRate = cond1_eeg.sampRate;
        
        % get n channels (to save later with TF files)
        %%% nElects = size(eeg.chanLabels,2);
        
        eegs = eeg.data(:,:,:); % get scalp EEG (drop EOG electrodes)
        %artInd = eeg.arf.artIndCleaned.'; % grab artifact rejection index
        tois = ismember(eeg.preTime:1000/Fs:eeg.postTime,em.time); nTimes = length(tois); % index time points for analysis.
        
        % %     %% TOM EDIT TO PROCESS WHOLE EPOCH (WTF DATA EPOCHED FROM -.5 to 2)
        % %     tois = ones(1,size(eegs,3)); nTimes = length(tois);
        
        % Remove rejected trials
        %eegs = eegs(~artInd,:,:);
        %posBin = posBin(~artInd);
        
        %Jordan Added
        %trial_indx = trial_indx(~artInd);
        trial_indx = trial_indx';
        
        em.nTrials = length(posBin); 
        nTrials = em.nTrials; % # of good trials
        
        %----------------------------------------------------------------------
        
        % Preallocate Matrices
        tf_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nChans); tf_total = tf_evoked;
        C2_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nBins,nChans); C2_total = C2_evoked;
        em.blocks = nan(nTrials,nIter);  % create em.block to save block assignments
        
        %Jordan added
        tf_evoked_cond1 = []; %nan(nFreqs,nIter,nSamps,nBlocks/2,nChans); 
        tf_evoked_cond2 = []; 
        tf_evoked_cond3 = [];
        tf_evoked_cond4 = [];
        tf_total_cond1 = []; 
        tf_total_cond2 = []
        tf_total_cond3 = [];
        tf_total_cond4 = [];
        
        C2_evoked_cond1 = [];%nan(nFreqs,nIter,nSamps,nBlocks/2,nBins,nChans);
        C2_total_cond1 = C2_evoked_cond1;
        C2_evoked_low = C2_evoked_cond1;
        C2_total_low = C2_evoked_low;
        
        C2_total_cond1Shifted = C2_total_cond1;
        C2_total_lowShifted = C2_total_low;
        %Jordan added for looking at central vs lateralized IEM responses
        horzMerid_locs = [1,8,4,5];
        vertMerid_locs = [2,3,6,7];
                        
        tfTotal_periphcond1 = []; tf_total_centralcond1 = tfTotal_periphcond1;
        tfTotal_periphcond2 = []; tf_total_centralcond2 = tfTotal_periphcond2;
        
        % TOM ADDED (needed to convert to EEGLAB format to use new eeglab filter)
        EEG.data = permute(eegs,[2,3,1]); % converts to chans x times x trials
        EEG.srate = Fs;
        EEG.trials = size(EEG.data,3);
        EEG.nbchan = size(EEG.data,1);
        EEG.pnts = size(EEG.data,2);
        
        % Loop through each frequency
        for f = 1:nFreqs
            tic % start timing frequency loop
            fprintf('Frequency %d out of %d\n', f, nFreqs)
            
            %% get no. of electrodes
            nElectrodes = size(eeg.data,2);
            em.nElectrodes = nElectrodes;
            disp(['nElecrodes changed to :' num2str(nElectrodes)])
            
            %% BUTTERWORTH FILTER
            filterorder = 3;
            type = 'bandpass';
            [z1,p1] = butter(filterorder, [freqs(f,1), freqs(f,2)]./(EEG.srate/2),type);
            %freqz(z1,p1,[],250)
            data = double(EEG.data);
            tempEEG = NaN(size(data,1),EEG.pnts,size(data,3));
            for x = 1:size(data,1) % loop through chans
                for y = 1:size(data,3) % loop through trials
                    dataFilt1 = filtfilt(z1,p1,data(x,:,y)); % was filtfilt
                    tempEEG(x,:,y) = dataFilt1; % tymp = chans x times x trials
                end
            end
            
            eegBand = [];
            eegBand = tempEEG;
            
            %% apply hilbert to each channel and epoch in turn (this should be correct)
            eegs = [];
            for j=1:size(tempEEG,1) % chan loop
                for i=1:size(tempEEG,3) % trial loop
                    eegs(i,j,:) = hilbert(squeeze(tempEEG(j,:,i)));
                end
            end
            
            % eegs is trials x elects x times
            fdata_evoked = eegs;
            fdata_total = abs(eegs).^2;
            
            % Loop through each iteration
            for iter = 1:nIter
                
                disp(['nIteration :' num2str(iter) ' of ' num2str(nIter)])
                %--------------------------------------------------------------------------
                % Assign trials to blocks (such that trials per position are
                % equated within blocks)
                %--------------------------------------------------------------------------
                
                % preallocate arrays
                blocks = nan(size(posBin));
                shuffBlocks = nan(size(posBin));
               
                % count number of trials within each position bin
                %clear binCnt
                
                %tom modified to work with WTF
%                 binCnt = [];
%                 for bin = 1:nBins
%                     binCnt(bin) = sum(posBin == bin);
%                 end
                
                
%                 binCnt = [];
%                 binCnt1 = [];
%                 binCnt2 = [];
%                 binCnt3 = [];
%                 binCnt4 = [];
%                 
%                 for bin = 1:nBins
%                     binCnt1(bin) = sum(posBin(trial_indx==1) == bin);
%                     binCnt2(bin) = sum(posBin(trial_indx==2) == bin);
%                     binCnt3(bin) = sum(posBin(trial_indx==3) == bin);
%                     binCnt4(bin) = sum(posBin(trial_indx==4) == bin);
%                 end
%                 
%                 binCnt = [binCnt1,binCnt2,binCnt3,binCnt4];
                
                % choose method of determining the location bin with the min
                % number of trials (1= per condition, 2= across conditions)
%                 if findMinType==1
%                     minCnt = min(binCnt); % # of trials for position bin with fewest trials
%                 elseif findMinType==2
%                     %if minCnt is based across all four conditions.  We get
%                     %the values from "Minimum_Location_Bin_MAT.mat" created
%                     %earlier
%                     if ismember(sn,1)
%                         minCnt=minLocsData.minLocBinAllConds(1);
%                     elseif ismember(sn,2)
%                         minCnt=minLocsData.minLocBinAllConds(2);
%                     elseif ismember(sn,3)
%                         minCnt=minLocsData.minLocBinAllConds(3);
%                     elseif ismember(sn,4)
%                         minCnt=minLocsData.minLocBinAllConds(4);
%                     elseif ismember(sn,5)
%                         minCnt=minLocsData.minLocBinAllConds(5);
%                     elseif ismember(sn,6)
%                         minCnt=minLocsData.minLocBinAllConds(6);
%                     elseif ismember(sn,7)
%                         minCnt=minLocsData.minLocBinAllConds(7);
%                     elseif ismember(sn,8)
%                         minCnt=minLocsData.minLocBinAllConds(8);
%                     elseif ismember(sn,10)
%                         minCnt=minLocsData.minLocBinAllConds(9);
%                     elseif ismember(sn,11)
%                         minCnt=minLocsData.minLocBinAllConds(10);
%                     elseif ismember(sn,12)
%                         minCnt=minLocsData.minLocBinAllConds(11);
%                     elseif ismember(sn,13)
%                         minCnt=minLocsData.minLocBinAllConds(12);
%                     elseif ismember(sn,14)
%                         minCnt=minLocsData.minLocBinAllConds(13);
%                     elseif ismember(sn,15)
%                         minCnt=minLocsData.minLocBinAllConds(14);
%                     elseif ismember(sn,16)
%                         minCnt=minLocsData.minLocBinAllConds(15);
%                     elseif ismember(sn,17)
%                         minCnt=minLocsData.minLocBinAllConds(16);
%                     elseif ismember(sn,18)
%                         minCnt=minLocsData.minLocBinAllConds(17);
%                     elseif ismember(sn,19)
%                         minCnt=minLocsData.minLocBinAllConds(18);
%                     elseif ismember(sn,20)
%                         minCnt=minLocsData.minLocBinAllConds(19);
%                     elseif ismember(sn,21)
%                         minCnt=minLocsData.minLocBinAllConds(20);
%                     elseif ismember(sn,22)
%                         minCnt=minLocsData.minLocBinAllConds(21);
%                     elseif ismember(sn,23)
%                         minCnt=minLocsData.minLocBinAllConds(22);
%                     elseif ismember(sn,24)
%                         minCnt=minLocsData.minLocBinAllConds(23);
%                     elseif ismember(sn,25)
%                         minCnt=minLocsData.minLocBinAllConds(24);
%                     elseif ismember(sn,26)
%                         minCnt=minLocsData.minLocBinAllConds(25);
%                     elseif ismember(sn,27)
%                         minCnt=minLocsData.minLocBinAllConds(26);
%                     elseif ismember(sn,28)
%                         minCnt=minLocsData.minLocBinAllConds(27);
%                     elseif ismember(sn,29)
%                         minCnt=minLocsData.minLocBinAllConds(28);
%                     elseif ismember(sn,30)
%                         minCnt=minLocsData.minLocBinAllConds(29);
%                     elseif ismember(sn,31)
%                         minCnt=minLocsData.minLocBinAllConds(30);
%                     elseif ismember(sn,32)
%                         minCnt=minLocsData.minLocBinAllConds(31);
%                     elseif ismember(sn,33)
%                         minCnt=minLocsData.minLocBinAllConds(32);
%                     elseif ismember(sn,34)
%                         minCnt=minLocsData.minLocBinAllConds(33);
%                     elseif ismember(sn,35)
%                         minCnt=minLocsData.minLocBinAllConds(34);
%                     end
%                 end
%                 
%                 % reduce minCnt by 1 trial to ensure that ALL location bins
%                 % have some degree of trial randomization per IEM iteration
%                 % (without this, the loc bin with the minimum no. of trials
%                 % will not be randomized at all).
%                 minCnt = minCnt-1;
                
                nPerBin = floor(minCnt/nBlocks); % max # of trials such that the # of trials for each bin can be equated within each block
                
                % shuffle trials
                shuffInd = randperm(nTrials)'; % create shuffle index
                shuffBin = posBin(shuffInd); % shuffle trial order
                
                shuffTrial_indx = trial_indx(shuffInd);
                
                % take the 1st nPerBin x nBlocks trials for each position bin.
                for bin = 1:nBins
                    idx = find(shuffBin == bin); % get index for trials belonging to the current bin
                    idx = idx(1:nPerBin*nBlocks); % drop excess trials
                    %x = repmat(1:nBlocks',nPerBin,1); shuffBlocks(idx) = x; % assign randomly order trials to blocks
                    

                    % Jordan Added
                    idx_cond1 = find(shuffBin == bin & shuffTrial_indx == 1);
                    idx_cond1 = idx_cond1(1:nPerBin*nBlocks);
                                        
                    idx_cond2 = find(shuffBin == bin & shuffTrial_indx == 2);
                    idx_cond2 = idx_cond2(1:nPerBin*nBlocks);
                    
                    idx_cond3 = find(shuffBin == bin & shuffTrial_indx == 3);
                    idx_cond3 = idx_cond3(1:nPerBin*nBlocks);
                    
                    idx_cond4 = find(shuffBin == bin & shuffTrial_indx == 4);
                    idx_cond4 = idx_cond4(1:nPerBin*nBlocks);
                    
                    fixed_TrainIdx = [idx_cond1;idx_cond2;idx_cond3;idx_cond4];
                    
                    % assign cond1 trials to first two blocks and exercise
                    % trials to second two blocks
                    if rem(length(idx_cond1),2)
                        first_blockN = round(length(idx_cond1)/2);
                        second_blockN = length(idx_cond1) - first_blockN;
                        x_cond1 = [ones(1,first_blockN),repmat(2,1,second_blockN)];
                    else
                        x_cond1 = repmat(1:2,1,length(idx_cond1)/2);
                    end
                    
                    if rem(length(idx_cond2),2)
                        first_blockN = round(length(idx_cond2)/2);
                        second_blockN = length(idx_cond2) - first_blockN;
                        x_cond2 = [repmat(3,1,first_blockN),repmat(4,1,second_blockN)];
                    else
                        x_cond2 = repmat(3:4,1,length(idx_cond2)/2);
                    end
                    
                    if rem(length(idx_cond3),2)
                        first_blockN = round(length(idx_cond3)/2);
                        second_blockN = length(idx_cond3) - first_blockN;
                        x_cond3 = [repmat(5,1,first_blockN),repmat(6,1,second_blockN)];
                    else
                        x_cond3 = repmat(5:6,1,length(idx_cond3)/2);
                    end
                    
                    if rem(length(idx_cond4),2)
                        first_blockN = round(length(idx_cond4)/2);
                        second_blockN = length(idx_cond4) - first_blockN;
                        x_cond4 = [repmat(7,1,first_blockN),repmat(8,1,second_blockN)];
                    else
                        x_cond4 = repmat(7:8,1,length(idx_cond4)/2);
                    end
                    
                    
                    
                    fixed_x = [x_cond1, x_cond2, x_cond3, x_cond4];
                    shuffBlocks(fixed_TrainIdx) = fixed_x;
                    
                end
                
                % unshuffle block assignment
                blocks(shuffInd) = shuffBlocks;
                
                %Jordan Added 
                
                %%%HERE - CHECK THIS IS CORRECT???%%%
                
                block_conMat = [blocks, shuffTrial_indx];
                blockCon_trialNum = [];
                for iBlock = 1:nBlocks
                    for iCon = 1:4
                        blockCon_trialNum(iBlock, iCon) = sum(block_conMat(:,1) == iBlock & block_conMat(:,2) == iCon);
                    end
                end
                
                % save block assignment
                em.blocks(:,iter) = blocks; % block assignment
                em.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
                
                %Jordan Added
                em.trial_Indx(:,iter) = trial_indx; %save condition assignment for each trial
             
                %-------------------------------------------------------------------------
                
                % Average data for each position bin across blocks
                posBins = 1:nBins;
                blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
                blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
                
                blockDat_evoked_cond1 = [];%nan(nBins*nBlocks/2,nElectrodes,nSamps);
                blockDat_total_cond1 = nan(size(blockDat_evoked_cond1));
                
                blockDat_evoked_cond2 = nan(size(blockDat_evoked_cond1));
                blockDat_total_cond2 = nan(size(blockDat_evoked_cond1));
                
                blockDat_evoked_cond3 = nan(size(blockDat_evoked_cond1));
                blockDat_total_cond3 = nan(size(blockDat_evoked_cond1));
                
                blockDat_evoked_cond4 = nan(size(blockDat_evoked_cond1));
                blockDat_total_cond4 = nan(size(blockDat_evoked_cond1));
                
                labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
                blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
                c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
                bCnt = 1; 
                
                % Jordan added
                cond1_labels = [];%nan(nBins*nBlocks/2,1);
                cond1_blockNum = nan(size(cond1_labels));
                cond1_c =[];%nan(nBins*nBlocks/2,nChans);
                
                cond2_labels = nan(size(cond1_labels));
                cond2_blockNum = nan(size(cond2_labels));
                cond2_c = nan(size(cond1_c));
                
                cond3_labels = nan(size(cond1_labels));
                cond3_blockNum = nan(size(cond2_labels));
                cond3_c = nan(size(cond1_c));
                
                cond4_labels = nan(size(cond1_labels));
                cond4_blockNum = nan(size(cond2_labels));
                cond4_c = nan(size(cond1_c));
                
                
                
                cond1_bCnt = 1;
                cond2_bCnt = 1;
                cond3_bCnt = 1;
                cond4_bCnt = 1;
                
                
                
                %%%%% HERE %%%%%
                
                
                for ii = 1:nBins
                    for iii = 1:nBlocks
                        blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                        blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1));
                        
                        
                        % Seperate condition data
                        if (iii == 1 || iii == 2) % cond1 condition
                            blockDat_evoked_cond1(cond1_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                            blockDat_total_cond1(cond1_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                            
                            
                            % Con labels for testing
                            cond1_labels(cond1_bCnt) = ii;
                            cond1_blockNum(cond1_bCnt) = iii;
                            cond1_c(cond1_bCnt,:) = basisSet(ii,:);
                            
                            cond1_bCnt = cond1_bCnt + 1;
                            
                            
                        elseif (iii == 3 || iii == 4) % Exercise condition
                            blockDat_evoked_cond2(cond2_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                            blockDat_total_cond2(cond2_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                        
                            cond2_labels(cond2_bCnt) = ii;
                            cond2_blockNum(cond2_bCnt) = iii;
                            cond2_c(cond2_bCnt,:) = basisSet(ii,:);
                            
                            cond2_bCnt = cond2_bCnt + 1;
                            
                        elseif (iii == 5 || iii == 6) % Exercise condition
                            blockDat_evoked_cond3(cond3_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                            blockDat_total_cond3(cond3_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                            
                            cond3_labels(cond3_bCnt) = ii;
                            cond3_blockNum(cond3_bCnt) = iii;
                            cond3_c(cond3_bCnt,:) = basisSet(ii,:);
                            
                            cond3_bCnt = cond3_bCnt + 1;
                            
                        elseif (iii == 7 || iii == 8) % Exercise condition
                            blockDat_evoked_cond4(cond4_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                            blockDat_total_cond4(cond4_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                            
                            cond4_labels(cond4_bCnt) = ii;
                            cond4_blockNum(cond4_bCnt) = iii;
                            cond4_c(cond4_bCnt,:) = basisSet(ii,:);
                            
                            cond4_bCnt = cond4_bCnt + 1;
                            
                            
                            
                        end
                        
                        
                        labels(bCnt) = ii; %provide labels for training and testing
                        blockNum(bCnt) = iii;
                        c(bCnt,:) = basisSet(ii,:);
                        bCnt = bCnt+1;
                        
                    end
                end
                
                for t = 1:nSamps
                    
                    % grab data for timepoint t
                    toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of intecond1
                    de = squeeze(mean(blockDat_evoked(:,:,toi),3)); % evoked data
                    dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
                    
                    %Jordan added
                    de_cond1 = squeeze(mean(blockDat_evoked_cond1(:,:,toi),3));
                    dt_cond1 = squeeze(mean(blockDat_total_cond1(:,:,toi),3));
                    
                    de_cond2 = squeeze(mean(blockDat_evoked_cond2(:,:,toi),3));
                    dt_cond2 = squeeze(mean(blockDat_total_cond2(:,:,toi),3));
                    
                    de_cond3 = squeeze(mean(blockDat_evoked_cond3(:,:,toi),3));
                    dt_cond3 = squeeze(mean(blockDat_total_cond3(:,:,toi),3));
                    
                    de_cond4 = squeeze(mean(blockDat_evoked_cond4(:,:,toi),3));
                    dt_cond4 = squeeze(mean(blockDat_total_cond4(:,:,toi),3));
                    
                    
                    % Do forward model
                    cond1Evoked_i = 1; cond2Evoked_i = 1; cond1Total_i = 1; cond2Total_i = 1;
                    cond3Evoked_i = 1; cond4Evoked_i = 1; cond3Total_i = 1; cond4Total_i = 1;

                    for i=1:nBlocks % loop through blocks, holding each out as the test set
                        
                        trnl = labels(blockNum~=i); % training labels
                        tstl = labels(blockNum==i); % test labels
                        
                        
                        cond1_tstl = cond1_labels(cond1_blockNum == i);
                        cond2_tstl = cond2_labels(cond2_blockNum == i);
                        cond3_tstl = cond3_labels(cond3_blockNum == i);
                        cond4_tstl = cond4_labels(cond4_blockNum == i);
                        
%                         %-----------------------------------------------------%
%                         % Analysis on Evoked Power                            %
%                         %-----------------------------------------------------%
%                         B1 = de(blockNum~=i,:);    % training data
%                         B2 = de(blockNum==i,:);    % test data
%                         C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
%                         W = C1\B1;          % estimate weight matrix
%                         C2 = (W'\B2')';     % estimate channel responses
%                         
%                         C2_evoked(f,iter,t,i,:,:) = C2; % save the unshifted channel responses
%                         
%                         % Save weight matrix (Mary MacLean added)
%                         WeightsTotal(:,:,i) = W;
%                         
%                         % shift eegs to common center
%                         n2shift = ceil(size(C2,2)/2);
%                         for ii=1:size(C2,1)
%                             [~, shiftInd] = min(abs(posBins-tstl(ii)));
%                             C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
%                         end
%                         
%                         tf_evoked(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
%                         
%                         % Jordan - after training on both, test on conditions
%                         % seperately
%                         if i == 1 || i == 2
%                             B2_cond1 = de_cond1(cond1_blockNum==i,:);
%                             C2_cond1 = (W'\B2_cond1')'; % test
%                             C2_evoked_cond1(f,iter,t,cond1Evoked_i,:,:) = C2_cond1;
%                             
%                             %shift
%                             n2shift = ceil(size(C2_cond1,2)/2);
%                             for ii=1:size(C2_cond1,1)
%                                 [~, shiftInd] = min(abs(posBins-cond1_tstl(ii)));
%                                 C2_cond1(ii,:) = wshift('1D', C2_cond1(ii,:), shiftInd-n2shift-1);
%                             end
%                             
%                             tf_evoked_cond1(f,iter,t,cond1Evoked_i,:) = mean(C2_cond1,1);
%                             
%                             cond1Evoked_i = cond1Evoked_i + 1;
%                         elseif i == 3 || i == 4
%                             B2_low = de_low(cond2_blockNum==i,:);
%                             C2_low = (W'\B2_low')'; % test
%                             C2_evoked_low(f,iter,t,cond2Evoked_i,:,:) = C2_low;
%                             
%                             n2shift = ceil(size(C2_low,2)/2);
%                             for ii=1:size(C2_low,1)
%                                 [~, shiftInd] = min(abs(posBins-cond2_tstl(ii)));
%                                 C2_low(ii,:) = wshift('1D', C2_low(ii,:), shiftInd-n2shift-1);
%                             end
%                             
%                             
%                             tf_evoked_low(f,iter,t,cond2Evoked_i,:) = mean(C2_low,1);
%                             cond2Evoked_i = cond2Evoked_i + 1;
%                         end
                        
                        
                        %-----------------------------------------------------%
                        % Analysis on Total Power                             %
                        %-----------------------------------------------------%
                        B1 = dt(blockNum~=i,:);    % training data
                        B2 = dt(blockNum==i,:);    % test data
                        C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
                        W = C1\B1;          % estimate weight matrix
                        C2 = (W'\B2')';     % estimate channel responses
                        
                        C2_total(f,iter,t,i,:,:) = C2;
                        
                        % shift eegs to common center
                        n2shift = ceil(size(C2,2)/2);
                        for ii=1:size(C2,1)
                            [~, shiftInd] = min(abs(posBins-tstl(ii)));
                            C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
                        end
                        
                        tf_total(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
                        
                        % Jordan - after training on both, test on conditions
                        % seperately
                        
                        if i == 1 || i == 2
                            B2_cond1 = dt_cond1(cond1_blockNum==i,:);
                            C2_cond1 = (W'\B2_cond1')'; % test
                            C2_total_cond1(f,iter,t,cond1Total_i,:,:) = C2_cond1;
                            
                            %shift conditions
                            n2shift = ceil(size(C2_cond1,2)/2);
                            for ii=1:size(C2_cond1,1)
                                [~, shiftInd] = min(abs(posBins-cond1_tstl(ii)));
                                C2_cond1(ii,:) = wshift('1D', C2_cond1(ii,:), shiftInd-n2shift-1);
                            end
                            
                            tf_total_cond1(f,iter,t,cond1Total_i,:) = mean(C2_cond1,1);
                            
                            C2_total_cond1Shifted(f,iter,t,cond1Total_i,:,:) = C2_cond1;
                            % Jordan added for JOCN reviewer to look at
                            % lateralized channel responses
                            
                            cond1Total_i = cond1Total_i + 1;
                            
                        elseif i == 3 || i == 4
                            
                            B2_cond2 = dt_cond2(cond2_blockNum==i,:);
                            C2_cond2 = (W'\B2_cond2')'; % test
                            C2_total_cond2(f,iter,t,cond2Total_i,:,:) = C2_cond2;
                            
                            
                            n2shift = ceil(size(C2_cond2,2)/2);
                            for ii=1:size(C2_cond2,1)
                                [~, shiftInd] = min(abs(posBins-cond2_tstl(ii)));
                                C2_cond2(ii,:) = wshift('1D', C2_cond2(ii,:), shiftInd-n2shift-1);
                            end
                            
                            
                            tf_total_cond2 (f,iter,t,cond2Total_i,:)  = mean(C2_cond2,1);
                            
                            C2_total_cond2Shifted(f,iter,t,cond2Total_i,:,:) = C2_cond2;
                            
                            cond2Total_i = cond2Total_i + 1;
                            
                            
                         elseif i == 5 || i == 6
                            
                            B2_cond3 = dt_cond3(cond3_blockNum==i,:);
                            C2_cond3 = (W'\B2_cond3')'; % test
                            C2_total_cond3(f,iter,t,cond3Total_i,:,:) = C2_cond3;
                            
                            
                            n2shift = ceil(size(C2_cond3,2)/2);
                            for ii=1:size(C2_cond3,1)
                                [~, shiftInd] = min(abs(posBins-cond3_tstl(ii)));
                                C2_cond3(ii,:) = wshift('1D', C2_cond3(ii,:), shiftInd-n2shift-1);
                            end
                            
                            
                            tf_total_cond3 (f,iter,t,cond3Total_i,:)  = mean(C2_cond3,1);
                            
                            C2_total_cond3Shifted(f,iter,t,cond3Total_i,:,:) = C2_cond3;
                            
                            cond3Total_i = cond3Total_i + 1;
                            
                          elseif i == 7 || i == 8
                            
                            B2_cond4 = dt_cond4(cond4_blockNum==i,:);
                            C2_cond4 = (W'\B2_cond4')'; % test
                            C2_total_cond4(f,iter,t,cond4Total_i,:,:) = C2_cond4;
                            
                            
                            n2shift = ceil(size(C2_cond4,2)/2);
                            for ii=1:size(C2_cond4,1)
                                [~, shiftInd] = min(abs(posBins-cond4_tstl(ii)));
                                C2_cond4(ii,:) = wshift('1D', C2_cond4(ii,:), shiftInd-n2shift-1);
                            end
                            
                            
                            tf_total_cond4 (f,iter,t,cond4Total_i,:)  = mean(C2_cond4,1);
                            
                            C2_total_cond4Shifted(f,iter,t,cond4Total_i,:,:) = C2_cond4;
                            
                            cond4Total_i = cond4Total_i + 1;
                            
                            
                            
                        end
                       
                        %-----------------------------------------------------%
                        
                    end
                    % Average weights across blocks add to matrix (Mary
                    % MacLean, added)
                    %allWeightsTotal(:,:,:,t,f,iter) = WeightsTotal;
                end
            end
            toc % stop timing the frequency loop
        end
        
        %% TOM ADDED
        % average over 1000 ITS + 3 BLOCKS to reduce size of saved file!
        %tf_evoked = squeeze(mean(mean(tf_evoked,2),4)); %tfs are nBins x ms 
        tf_total = squeeze(mean(mean(tf_total,2),4));
        
        %Jordan Added
        %tf_evoked_cond1 = squeeze(mean(mean(tf_evoked_cond1,2),4));
        %tf_evoked_cond2 =  squeeze(mean(mean(tf_evoked_cond2,2),4));
        
        tf_total_cond1 = squeeze(mean(mean(tf_total_cond1,2),4));
        tf_total_cond2 =  squeeze(mean(mean(tf_total_cond2,2),4));
        tf_total_cond3 = squeeze(mean(mean(tf_total_cond3,2),4));
        tf_total_cond4 =  squeeze(mean(mean(tf_total_cond4,2),4));
        
        
%         %Jordan added for JOCN Reviewer
%         tfTotal_periphcond1 = squeeze(mean(mean(tfTotal_periphcond1,2),4));
%         tf_total_centralcond1 = squeeze(mean(mean(tf_total_centralcond1,2),4));
%         
%         tfTotal_periphcond2 = squeeze(mean(mean(tfTotal_periphcond2,2),4));
%         tfTotal_centralcond2 = squeeze(mean(mean(tfTotal_centralcond2,2),4));
        
        cd([root '/Analysis_Scripts'])
        % save data
        if processForIEM==1
            %fName = [root, '/Train4Conds/',sprintf('sj%02d_Train4Conds_WTF_Test', subs(s)),name];
            fName = [root, '/IEM_Results_TT_Within_Fixed/',sprintf('sj%02d_fixed_IEM.mat',subs(s))];
            
            %em.C2.evoked = C2_evoked;
            em.C2.total = C2_total;
            %em.tfs.evoked = tf_evoked;
            em.tfs.total = tf_total;
            %em.tfs.totalW = allWeightsTotal;
            em.nBlocks = nBlocks;
            %save(fName,'em','-v7.3');
            %em.C2_cond1.evoked = C2_evoked_cond1;
            
            %em.C2_cond1.total = C2_total_cond1;
            %em.C2_low.evoked = C2_evoked_low;
            %em.C2_cond2.total = C2_total_cond2;
            %em.C2_cond3.total = C2_total_cond3;
            %em.C2_cond4.total = C2_total_cond4;

            
            %em.tfs_cond1.evoked = tf_evoked_cond1;
            em.tfs_cond1.total = tf_total_cond1;
            %em.tfs_low.evoked = tf_evoked_low;
            em.tfs_cond2.total = tf_total_cond2;
            em.tfs_cond3.total = tf_total_cond3;
            em.tfs_cond4.total = tf_total_cond4;
            
%             em.splitTFs.cond1.total.peripheral = tfTotal_periphcond1;
%             em.splitTFs.cond1.total.central = tfTotal_centralcond1;
%             em.splitTFs.cond2.total.peripheral = tfTotal_periphcond2;
%             em.splitTFs.cond2.total.central = tfTotal_centralcond2;
            
            %em.C2_cond1.shifted.total = C2_total_cond1Shifted;
            %em.C2_low.shitfed.total = C2_total_lowShifted;
            
            %parsave(fName,em,minCnt,nElectrodes);
            save(fName,'em','minCnt','nElectrodes');
        else
            % save raw hilbert files
            fName = [hRoot,sprintf('sj%02d_TrainBoth_changeDect_accTrials', subs(s)),hilbName];
            eegInfo.chanLabels = eeg.chanLabels;
            eegInfo.preTime = eeg.preTime;
            eegInfo.postTime = eeg.postTime;
            eegInfo.sampRate = eeg.sampRate;
            eegInfo.posBin = posBin;
            parsave(fName, eegs, eegInfo, eegBand)
        end
        
    end
    
    %matlabpool close
    
end

%sendEmailToMe('SPATIAL IEM SCRIPT FINISHED PROCESSING!!')

clear all
close all
end