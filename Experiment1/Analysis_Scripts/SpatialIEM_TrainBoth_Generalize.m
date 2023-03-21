function SpatialIEM_TrainBoth_Generalize (subs)
%==========================================================================
%{
Purpose: run IEM with training/testing on all time samples within each of
the four conditions

Original Author:
Joshua J. Foster
joshua.james.foster@gmail.com
University of Chicago
August 12, 2015

Modified by Tom Bullock
UCSB Attention Lab
%}
%==========================================================================

%clear all
close all

%subs = [1:8,10];


nSubs = length(subs);

%% process for IEM (1=IEM)
processForIEM=1;

% setup directories
root = '/home/garrett/WTF_Bike/';
dRoot = [root, 'Data/'];

if processForIEM==1
    eRoot = [root,'EEG/'];
else
    eRoot = [root,'HILBERT/'];
end

hRoot = [root,'HILBERT/'];

% choose method for finding the min number of trials per location bin
% (1=per-condition, 2=across conditions)
findMinType = 2;

% load matrix of min trials per loc bin per condition
minLocsData = load([root 'Analysis_Scripts/Minimum_Location_Bin_Mat_AccTrials.mat']);



name = '_SpatialTF_ALPHA_Gen.mat'; % name of files to be saved for IEM

thisFreq = {'Alpha'};
freqBandpass = [8 12];



% Loop through participants (***THIS WAS DISABLED)
%  matlabpool open 72

em = [];
EEG = [];
eegInfo = [];

allB1=[];
allB2=[];
allC1=[];

% parameters to set
em.nChans = 8; % # of channels
em.nBins = em.nChans; % # of stimulus bins
em.nIter = 10; % # of iterations %%% WAS SET TO 10!!!
em.nBlocks = 4; % # of blocks for cross-validation
em.frequencies = freqBandpass; % frequency bands to analyze
em.bands = thisFreq;
em.window = 4;
em.Fs = 250; % WAS 250
em.nElectrodes = 999; % get later
em.time = -.5*1000:1000/em.Fs:1.9961*1000; %   -500:4:2000; % time points of interest

% determine downsample rate for generalizations
em.genSamps = length(em.time)/25; %divide 625 samples into 25 samples of 25 (100ms per sample)

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
genSamps = em.genSamps; %generalization samples

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

fprintf('Subject:\t%d\n',subs)

% Grab data------------------------------------------------------------

% Get position bin index from behavior file
restfName = [dRoot, sprintf('sj%02d_exerCon01_changeDect_MixModel_wBias_accTrials.mat',subs)];
exfName = [dRoot, sprintf('sj%02d_exerCon02_changeDect_MixModel_wBias_accTrials.mat',subs)];

% load files
rest_tmp = []; rest_beh = [];
rest_tmp = load(restfName);
rest_beh = rest_tmp.beh;

ex_tmp = []; ex_beh = [];
ex_tmp = load(exfName);
ex_beh = ex_tmp.beh;

% Position info
beh.trial.posBin = [rest_beh.trial.posBin,ex_beh.trial.posBin]; %merge data

%Jordan Added: create index to keep track of rest and exercise position trials
trial_indx = [repmat(1,1,length(rest_beh.trial.posBin)), repmat(2,1,length(ex_beh.trial.posBin))];

em.posBin = beh.trial.posBin'; % add to fm structure so it's saved
posBin = em.posBin;

% Performance info
%em.error = beh.trial.err';
%errorData = [beh.trial.err];
%medianError = median(errorData);

% Get EEG data
restfName = [eRoot, sprintf('sj%02d_exerCon01_changeDect_EEG_accTrials.mat',subs)];
exfName = [eRoot, sprintf('sj%02d_exerCon02_changeDect_EEG_accTrials.mat',subs)];


% load file
rest_tmp = []; rest_eeg = [];
rest_tmp = load(restfName);
rest_eeg = rest_tmp.eeg;

ex_tmp = []; ex_eeg = [];
ex_tmp = load(exfName);
ex_eeg = ex_tmp.eeg;

%Jordan Added
eeg.data = cat(1,rest_eeg.data,ex_eeg.data); %merge eeg trials along first dimension (trials x channels x timepoints)
eeg.arf.artIndCleaned = [rest_eeg.arf.artIndCleaned, ex_eeg.arf.artIndCleaned];

eeg.preTime = rest_eeg.preTime;
eeg.postTime = rest_eeg.postTime;
eeg.chanLabels = rest_eeg.chanLabels;
eeg.sampRate = rest_eeg.sampRate;

eegs = eeg.data(:,:,:); % get scalp EEG (drop EOG electrodes)
artInd = eeg.arf.artIndCleaned.'; % grab artifact rejection index
tois = ismember(eeg.preTime:1000/Fs:eeg.postTime,em.time); nTimes = length(tois); % index time points for analysis.

% Remove rejected trials
eegs = eegs(~artInd,:,:);
posBin = posBin(~artInd);
%errorData = errorData(~artInd)';

%Jordan Added
trial_indx = trial_indx(~artInd);
trial_indx = trial_indx';

em.nTrials = length(posBin); 
nTrials = em.nTrials; % # of good trials


%----------------------------------------------------------------------
em.blocks = nan(nTrials,nIter);  % create em.block to save block assignments

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
    
    % make n-electrodes flexible to account for chan rej
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
    for x = 1:size(data,1)
        for y = 1:size(data,3)
            dataFilt1 = filtfilt(z1,p1,data(x,:,y)); % was filtfilt
            tempEEG(x,:,y) = dataFilt1; % tymp = chans x times x trials
        end
    end
    
    % apply hilbert to each channel and epoch in turn (this should be correct)
    eegs = [];
    for j=1:size(tempEEG,1) % chan loop
        for i=1:size(tempEEG,3) % trial loop
            eegs(i,j,:) = hilbert(squeeze(tempEEG(j,:,i)));
        end
    end
    
    % eegs is trials x elects x times
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
        binCnt = [];
        for bin = 1:nBins
            binCnt(bin) = sum(posBin == bin);
        end
        
        % choose method of determining the location bin with the min
        % number of trials (1= per condition, 2= across conditions)
        if findMinType==1
            minCnt = min(binCnt); % # of trials for position bin with fewest trials
        elseif findMinType==2
            %if minCnt is based across all four conditions.  We get
            %the values from "Minimum_Location_Bin_MAT.mat" created
            %earlier
            if ismember(subs,1)
                minCnt=minLocsData.minLocBinAllConds(1);
            elseif ismember(subs,2)
                minCnt=minLocsData.minLocBinAllConds(2);
            elseif ismember(subs,3)
                minCnt=minLocsData.minLocBinAllConds(3);
            elseif ismember(subs,4)
                minCnt=minLocsData.minLocBinAllConds(4);
            elseif ismember(subs,5)
                minCnt=minLocsData.minLocBinAllConds(5);
            elseif ismember(subs,6)
                minCnt=minLocsData.minLocBinAllConds(6);
            elseif ismember(subs,7)
                minCnt=minLocsData.minLocBinAllConds(7);
            elseif ismember(subs,8)
                minCnt=minLocsData.minLocBinAllConds(8);
            elseif ismember(subs,10)
                minCnt=minLocsData.minLocBinAllConds(9);
            elseif ismember(subs,11)
                minCnt=minLocsData.minLocBinAllConds(10);
            elseif ismember(subs,12)
                minCnt=minLocsData.minLocBinAllConds(11);
            elseif ismember(subs,13)
                minCnt=minLocsData.minLocBinAllConds(12);
            elseif ismember(subs,14)
                minCnt=minLocsData.minLocBinAllConds(13);
            elseif ismember(subs,15)
                minCnt=minLocsData.minLocBinAllConds(14);
            elseif ismember(subs,16)
                minCnt=minLocsData.minLocBinAllConds(15);
            elseif ismember(subs,17)
                minCnt=minLocsData.minLocBinAllConds(16);
            elseif ismember(subs,18)
                minCnt=minLocsData.minLocBinAllConds(17);
            elseif ismember(subs,19)
                minCnt=minLocsData.minLocBinAllConds(18);
            elseif ismember(subs,20)
                minCnt=minLocsData.minLocBinAllConds(19);
            elseif ismember(subs,21)
                minCnt=minLocsData.minLocBinAllConds(20);
            elseif ismember(subs,22)
                minCnt=minLocsData.minLocBinAllConds(21);
            elseif ismember(subs,23)
                minCnt=minLocsData.minLocBinAllConds(22);
            elseif ismember(subs,24)
                minCnt=minLocsData.minLocBinAllConds(23);
            elseif ismember(subs,25)
                minCnt=minLocsData.minLocBinAllConds(24);
            elseif ismember(subs,26)
                minCnt=minLocsData.minLocBinAllConds(25);
            elseif ismember(subs,27)
                minCnt=minLocsData.minLocBinAllConds(26);
            elseif ismember(subs,28)
                minCnt=minLocsData.minLocBinAllConds(27);
            elseif ismember(subs,29)
                minCnt=minLocsData.minLocBinAllConds(28);
            elseif ismember(subs,30)
                minCnt=minLocsData.minLocBinAllConds(29);
            elseif ismember(subs,31)
                minCnt=minLocsData.minLocBinAllConds(30);
            elseif ismember(subs,32)
                minCnt=minLocsData.minLocBinAllConds(31);
            elseif ismember(subs,33)
                minCnt=minLocsData.minLocBinAllConds(32);
            elseif ismember(subs,34)
                minCnt=minLocsData.minLocBinAllConds(33);
            elseif ismember(subs,35)
                minCnt=minLocsData.minLocBinAllConds(34);
                
            end
        end
        
        % reduce minCnt by 1 trial to ensure that ALL location bins
        % have some degree of trial randomization per IEM iteration
        % (without this, the loc bin with the minimum no. of trials
        % will not be randomized at all).
        minCnt = minCnt-1;
        
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
            idx_rest = find(shuffBin == bin & shuffTrial_indx == 1);
            idx_rest = idx_rest(1:nPerBin*nBlocks);
            
            idx_low = find(shuffBin == bin & shuffTrial_indx == 2);
            idx_low = idx_low(1:nPerBin*nBlocks);
            con_idx = shuffTrial_indx(idx);
            con_idx = con_idx(1:nPerBin*nBlocks)';
            
            
            fixed_TrainIdx = [idx_rest;idx_low];
            
            
            % assign rest trials to first two blocks and exercise
            % trials to second two blocks
            if rem(length(idx_rest),2)
                first_blockN = round(length(idx_rest)/2);
                second_blockN = length(idx_rest) - first_blockN;
                x_rest = [ones(1,first_blockN),repmat(2,1,second_blockN)];
            else
                x_rest = repmat(1:2,1,length(idx_rest)/2);
            end
            
            if rem(length(idx_low),2)
                first_blockN = round(length(idx_low)/2);
                second_blockN = length(idx_low) - first_blockN;
                x_low = [repmat(3,1,first_blockN),repmat(4,1,second_blockN)];
            else
                x_low = repmat(3:4,1,length(idx_low)/2);
            end
            
            fixed_x = [x_rest, x_low];
            shuffBlocks(fixed_TrainIdx) = fixed_x;
            
        end
        
        
        % unshuffle block assignment
        blocks(shuffInd,:) = shuffBlocks;
        
        %Jordan Added
        
        block_conMat = [blocks, shuffTrial_indx];
        blockCon_trialNum = [];
        for iBlock = 1:nBlocks
            for iCon = 1:2
                blockCon_trialNum(iBlock, iCon) = sum(block_conMat(:,1) == iBlock & block_conMat(:,2) == iCon);
            end
        end
        
        % save block assignment
        em.blocks(:,iter) = blocks; % block assignment
        em.nTrialsPerBlock = length(blocks(blocks == 1)); % # of trials per block
        
        %Jordan Added
        em.trial_Indx(:,iter) = trial_indx;
        
        %-------------------------------------------------------------------------
        
        % Average data for each position bin across blocks
        posBins = 1:nBins;
        blockDat_evoked = nan(nBins*nBlocks,nElectrodes,nSamps); % averaged evoked data
        blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
        
        blockDat_evoked_Rest = [];%nan(nBins*nBlocks/2,nElectrodes,nSamps);
        blockDat_total_Rest = nan(size(blockDat_evoked_Rest));
        
        blockDat_evoked_Low = nan(size(blockDat_evoked_Rest));
        blockDat_total_Low = nan(size(blockDat_evoked_Low));
        
        labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
        blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
        c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
        bCnt = 1;
        
        % Jordan added
        rest_labels = [];%nan(nBins*nBlocks/2,1);
        rest_blockNum = nan(size(rest_labels));
        rest_c =[];%nan(nBins*nBlocks/2,nChans);
        
        ex_labels = nan(size(rest_labels));
        ex_blockNum = nan(size(ex_labels));
        ex_c = nan(size(rest_c));
        
        rest_bCnt = 1;
        ex_bCnt = 1;
        for ii = 1:nBins
            for iii = 1:nBlocks
                
                blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1));
                
                
                % Seperate condition data
                if (iii == 1 || iii == 2) % Rest condition
                    
                    blockDat_total_Rest(rest_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                    
                    
                    % Con labels for testing
                    rest_labels(rest_bCnt) = ii;
                    rest_blockNum(rest_bCnt) = iii;
                    rest_c(rest_bCnt,:) = basisSet(ii,:);
                    
                    rest_bCnt = rest_bCnt + 1;
                elseif (iii == 3 || iii == 4) % Exercise condition
                    
                    blockDat_total_Low(ex_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1)));
                    
                    ex_labels(ex_bCnt) = ii;
                    ex_blockNum(ex_bCnt) = iii;
                    ex_c(ex_bCnt,:) = basisSet(ii,:);
                    
                    ex_bCnt = ex_bCnt + 1;
                end
                
                labels(bCnt) = ii; %provide labels for training and testing
                blockNum(bCnt) = iii;
                c(bCnt,:) = basisSet(ii,:);
                bCnt = bCnt+1;
            end
        end
        
        
        % training loop
        tf_total_rest=[];
        tf_total_ex = [];
        for trLoop=1:genSamps %divide 625 samples into 25 samples of 25 (100ms per sample)
            
            this_trainSamps = (genSamps*trLoop)-(genSamps-1):(genSamps*trLoop);% divide 640 samples into 40 samples of 16 (62.5ms per sample)
            
            trainData = squeeze(mean(blockDat_total(:,:,this_trainSamps),3));
            
            
            % testing loop
            for teLoop=1:genSamps
                
                this_testSamps = (genSamps*teLoop)-(genSamps-1):(genSamps*teLoop);
                testData = squeeze(mean(blockDat_total(:,:,this_testSamps),3));
                
                % Do forward model
                restTotal_i = 1; exTotal_i = 1;
                for i=1:nBlocks % loop through blocks, holding each out as the test set
                    
                    
                    trnl = labels(blockNum~=i); % training labels
                    
                    rest_tstl = rest_labels(rest_blockNum == i);
                    ex_tstl = ex_labels(ex_blockNum == i);
                    
                    %-----------------------------------------------------%
                    % Analysis on Total Power                             %
                    %-----------------------------------------------------%
                    
                    B1 = trainData(blockNum~=i,:);    % training data
                    C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
                    W = C1\B1;
                    
                    if i == 1 || i == 2
                        rest_testData = squeeze(mean(blockDat_total_Rest(:,:,this_testSamps),3));
                        B2_rest = rest_testData(rest_blockNum == i,:);
                        C2_rest = (W'\B2_rest')';
                        
                        %shift
                        n2shift = ceil(size(C2_rest,2)/2);
                        for ii=1:size(C2_rest,1)
                            [~, shiftInd] = min(abs(posBins-rest_tstl(ii)));
                            C2_rest(ii,:) = wshift('1D', C2_rest(ii,:), shiftInd-n2shift-1);
                        end
                        
                        tf_total_rest(iter,trLoop,teLoop,restTotal_i,:) = mean(C2_rest,1);
                        restTotal_i = restTotal_i + 1;
                        
                    elseif i == 3 || i == 4
                        ex_testData = squeeze(mean(blockDat_total_Low(:,:,this_testSamps),3));
                        B2_ex = ex_testData(ex_blockNum == i,:);
                        C2_ex = (W'\B2_ex')';
                        
                        %shift
                        n2shift = ceil(size(C2_ex,2)/2);
                        for ii=1:size(C2_ex,1)
                            [~, shiftInd] = min(abs(posBins-ex_tstl(ii)));
                            C2_ex(ii,:) = wshift('1D', C2_ex(ii,:), shiftInd-n2shift-1);
                        end
                        
                        tf_total_ex(iter,trLoop,teLoop,exTotal_i,:) = mean(C2_ex,1);
                        exTotal_i = exTotal_i + 1;
                    end
                end
            end
        end
    end
end

toc % stop timing the frequency loop

% average over iterations and blocks to reduce file

tf_total_rest = squeeze(mean(mean(tf_total_rest,1),4));
tf_total_ex = squeeze(mean(mean(tf_total_ex,1),4));

% Jordan added to look at raw data
em.raw_total = fdata_total;
% save data
if processForIEM==1 % IEM
    
    fName = [dRoot,'TrainBoth/', sprintf('sj%02d_TrainBoth_changeDect_acc',subs), name];
    
    %em.tfs.evoked = tf_evoked;
    em.tfs.total.rest = tf_total_rest;
    em.tfs.total.ex = tf_total_ex;
    em.nBlocks = nBlocks;
    em.c=c;
    %save(fName,'em','-v7.3');
    save(fName,'em','minCnt','nElectrodes','-v7.3');
else %HILBERT
    fName = [hRoot,num2str(subs),hilbName];
    eegInfo.chanLabels = eeg.chanLabels;
    eegInfo.preTime = eeg.preTime;
    eegInfo.postTime = eeg.postTime;
    eegInfo.sampRate = eeg.sampRate;
    eegInfo.posBin = posBin;
    parsave(fName, eegs, eegInfo)
end


close all
clear all


%sendEmailToMe('SPATIAL IEM SCRIPT FINISHED PROCESSING!!')
