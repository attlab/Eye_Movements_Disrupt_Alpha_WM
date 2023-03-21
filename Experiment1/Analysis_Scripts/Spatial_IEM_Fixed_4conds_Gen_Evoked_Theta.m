function Spatial_IEM_Fixed_4conds_Gen_Evoked_Theta(sjNum)

%{
Purpose: run spatial encoding fixed model 

Original Author:
Joshua J. Foster
joshua.james.foster@gmail.com
University of Chicago
August 12, 2015
Modified heavily by Jordan Garrett and Tom Bullock, 
UCSB Attention Lab
%}

% process for IEM (1) or Hilbert (2)
processForIEM=1;

% setup directories
root = '/home/waldrop/Desktop/WTF_EYE';
eegDir = [root '/' 'EEG_Prepro2_Avg_Baseline'];

% frequency band settings (may drop)
thisFreq = {'Theta'};
freqBandpass = [4 7];

em = [];
EEG = [];
eegInfo = [];

% parameters to set
em.nChans = 8; % # of hypothetical ori/loc channels;
em.nBins = em.nChans; % # of stimulus bins
em.nIter = 10; % # of iterations 
em.nBlocks = 8; % # of blocks for cross-validation
em.frequencies = freqBandpass; % frequency bands to analyze
em.bands = thisFreq;
em.window = 4;
em.Fs = 256;
em.nElectrodes = 999; % get later
em.time = -.5*1000:1000/em.Fs:1.9961*1000; % -500:4:2000; % time points of interest

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


% merge EEG/BEH data across sessions and split by condition (SCAM)
for iSession=1:2
    
    % load data
    load([eegDir '/' sprintf('sj%02d_se%02d_wtfEye.mat',sjNum,iSession)])
    
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

% parse behavioral data again
cond1_beh = [allConds(1).beh.stimLoc];
cond2_beh = [allConds(2).beh.stimLoc];
cond3_beh = [allConds(3).beh.stimLoc];
cond4_beh = [allConds(4).beh.stimLoc];

clear beh
% merge position data into a single vector
beh.trial.posBin = [cond1_beh,cond2_beh,cond3_beh,cond4_beh];

% creat index to keep track of conditions [RENAME]
trial_indx = [repmat(1,1,length(cond1_beh)),repmat(2,1,length(cond2_beh)),repmat(3,1,length(cond3_beh)),repmat(4,1,length(cond4_beh))];

em.posBin = beh.trial.posBin'; % add to fm structure so it's saved
clear posBin
posBin = em.posBin;

% organize eeg data
cond1_eeg = allConds(1).eegs; % chans x times x trials
cond2_eeg = allConds(2).eegs; % chans x times x trials
cond3_eeg = allConds(3).eegs; % chans x times x trials
cond4_eeg = allConds(4).eegs; % chans x times x trials

% permute to correct format (trials x chans x times)
cond1_eeg = permute(cond1_eeg,[3,1,2]);
cond2_eeg = permute(cond2_eeg,[3,1,2]);
cond3_eeg = permute(cond3_eeg,[3,1,2]);
cond4_eeg = permute(cond4_eeg,[3,1,2]);

% recreate eeg structure
eeg.preTime = eeg.preTime;
eeg.postTime = eeg.postTime;
eeg.chanLabels = eeg.chanLabels;
eeg.sampRate = eeg.sampRate;

eeg.data = [];
eeg.data = cat(1,cond1_eeg,cond2_eeg,cond3_eeg,cond4_eeg);

% isolate eeg data
eegs = eeg.data(:,:,:);

% set times of interest (i.e. timepoints that get run through the IEM)
tois = ismember(eeg.preTime:1000/Fs:eeg.postTime,em.time); 

% transpose trial_indx to a column vector
trial_indx = trial_indx';

em.nTrials = length(posBin);
nTrials = em.nTrials; % # of good trials

%----------------------------------------------------------------------

% Preallocate Matrices
tf_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nChans); tf_total = tf_evoked;
C2_evoked = nan(nFreqs,nIter,nSamps,nBlocks,nBins,nChans); C2_total = C2_evoked;
em.blocks = nan(nTrials,nIter);  % create em.block to save block assignments

% set up some vars (probably not necessary)
tf_total_cond1 = [];
tf_total_cond2 = [];
tf_total_cond3 = [];
tf_total_cond4 = [];

C2_total_cond1 = [];
C2_total_cond2 = [];
C2_total_cond3 = [];
C2_total_cond4 = [];

% convert eegs to appropriate format for Butterworth filter
EEG.data = permute(eegs,[2,3,1]); % converts to chans x times x trials
EEG.srate = Fs;
EEG.trials = size(EEG.data,3);
EEG.nbchan = size(EEG.data,1);
EEG.pnts = size(EEG.data,2);

% Loop through each frequency
for f = 1:nFreqs
    tic % start timing frequency loop
    fprintf('Frequency %d out of %d\n', f, nFreqs)
    
    % determine number of electrodes
    nElectrodes = size(eeg.data,2);
    em.nElectrodes = nElectrodes;
    disp(['nElecrodes changed to :' num2str(nElectrodes)])
    
    % run Butterworth filter
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
    
    % apply hilbert to each channel and epoch in turn 
    eegs = [];
    for j=1:size(tempEEG,1) % chan loop
        for i=1:size(tempEEG,3) % trial loop
            eegs(i,j,:) = hilbert(squeeze(tempEEG(j,:,i)));
        end
    end
    
    % eegs is trials x elects x times
    %fdata_total = abs(eegs).^2;
    fdata_total = eegs;
    
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
        
        % max # trials per bin so that blocks can be equated
        nPerBin = floor(minCnt/nBlocks); 
        
        % shuffle trials
        shuffInd = randperm(nTrials)'; % create shuffle index
        shuffBin = posBin(shuffInd); % shuffle trial order
        shuffTrial_indx = trial_indx(shuffInd);
        
        % take the 1st nPerBin x nBlocks trials for each position bin.
        for bin = 1:nBins
            idx = find(shuffBin == bin); % get index for trials belonging to the current bin
            idx = idx(1:nPerBin*nBlocks); % drop excess trials
        
            % invididual conditions
            idx_cond1 = find(shuffBin == bin & shuffTrial_indx == 1);
            idx_cond1 = idx_cond1(1:nPerBin*nBlocks);
            
            idx_cond2 = find(shuffBin == bin & shuffTrial_indx == 2);
            idx_cond2 = idx_cond2(1:nPerBin*nBlocks);
            
            idx_cond3 = find(shuffBin == bin & shuffTrial_indx == 3);
            idx_cond3 = idx_cond3(1:nPerBin*nBlocks);
            
            idx_cond4 = find(shuffBin == bin & shuffTrial_indx == 4);
            idx_cond4 = idx_cond4(1:nPerBin*nBlocks);
            
            fixed_TrainIdx = [idx_cond1;idx_cond2;idx_cond3;idx_cond4];
            
            % assign the conditions to pairs of blocks
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
        
        % set up blocks
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
        blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
        
        blockDat_total_cond1 = [];
        blockDat_total_cond2 = [];
        blockDat_total_cond3 = [];
        blockDat_total_cond4 = [];
        
        labels = nan(nBins*nBlocks,1); % bin labels for averaged data
        blockNum = nan(nBins*nBlocks,1); % block numbers for averaged data
        c = nan(nBins*nBlocks,nChans); % predicted channel responses for averaged data
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
        
        
        for ii = 1:nBins
            for iii = 1:nBlocks
                %blockDat_evoked(bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                blockDat_total(bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                
                
                % Seperate condition data
                if (iii == 1 || iii == 2) % cond1 condition
                    %blockDat_evoked_cond1(cond1_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    blockDat_total_cond1(cond1_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    
                    cond1_labels(cond1_bCnt) = ii;
                    cond1_blockNum(cond1_bCnt) = iii;
                    cond1_c(cond1_bCnt,:) = basisSet(ii,:);
                    cond1_bCnt = cond1_bCnt + 1;
                    
                elseif (iii == 3 || iii == 4) % Exercise condition
                    %blockDat_evoked_cond2(cond2_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    blockDat_total_cond2(cond2_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    
                    cond2_labels(cond2_bCnt) = ii;
                    cond2_blockNum(cond2_bCnt) = iii;
                    cond2_c(cond2_bCnt,:) = basisSet(ii,:);
                    cond2_bCnt = cond2_bCnt + 1;
                    
                elseif (iii == 5 || iii == 6) % Exercise condition
                    %blockDat_evoked_cond3(cond3_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    blockDat_total_cond3(cond3_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    
                    cond3_labels(cond3_bCnt) = ii;
                    cond3_blockNum(cond3_bCnt) = iii;
                    cond3_c(cond3_bCnt,:) = basisSet(ii,:);
                    
                    cond3_bCnt = cond3_bCnt + 1;
                    
                elseif (iii == 7 || iii == 8) % Exercise condition
                    %blockDat_evoked_cond4(cond4_bCnt,:,:) = abs(squeeze(mean(fdata_evoked(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    blockDat_total_cond4(cond4_bCnt,:,:) = abs(squeeze(mean(fdata_total(posBin==posBins(ii) & blocks==iii,:,tois),1))).^2;
                    
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
        
        % START GENERALIZATION
        
        for trLoop = 1:40 %nSamps
            
           % get times of interest (downsampled)
            this_trainSamps = (16*trLoop)-15:(16*trLoop);
            
            
            %trainData = squeeze(mean(blockDat_total(:,:,this_trainSamps),3));
            
            % TRAINING DATA MAT
            dt = squeeze(mean(blockDat_total(:,:,this_trainSamps),3));  % total data
            
            
            
            
            for teLoop = 1:40
                
                
                % get times of interest (downsampled)
                this_testSamps = (16*teLoop)-15:(16*teLoop);
                
                
                
                
                %trainData = squeeze(mean(blockDat_total(:,:,theseSamples),3));
                
                % grab data for timepoint t
                %toi = ismember(times,times(t)-em.window/2:times(t)+em.window/2); % time window of intecond1
                %            de = squeeze(mean(blockDat_evoked(:,:,toi),3)); % evoked data
                %dt = squeeze(mean(blockDat_total(:,:,toi),3));  % total data
                
                % THESE ARE THE TEST DATA MATS
                dt_cond1 = squeeze(mean(blockDat_total_cond1(:,:,this_testSamps),3));
                dt_cond2 = squeeze(mean(blockDat_total_cond2(:,:,this_testSamps),3));
                dt_cond3 = squeeze(mean(blockDat_total_cond3(:,:,this_testSamps),3));
                dt_cond4 = squeeze(mean(blockDat_total_cond4(:,:,this_testSamps),3));
                
                % Do forward model
                cond1Total_i = 1; cond2Total_i = 1;cond3Total_i = 1; cond4Total_i = 1;
                
                for i=1:nBlocks % loop through blocks, holding each out as the test set
                    
                    trnl = labels(blockNum~=i); % training labels
                    tstl = labels(blockNum==i); % test labels
                    
                    cond1_tstl = cond1_labels(cond1_blockNum == i);
                    cond2_tstl = cond2_labels(cond2_blockNum == i);
                    cond3_tstl = cond3_labels(cond3_blockNum == i);
                    cond4_tstl = cond4_labels(cond4_blockNum == i);
                    
                    %-----------------------------------------------------%
                    % Analysis on Total Power                             %
                    %-----------------------------------------------------%
                    B1 = dt(blockNum~=i,:);    % training data
                    B2 = dt(blockNum==i,:);    % test data
                    C1 = c(blockNum~=i,:);     % predicted channel outputs for training data
                    W = C1\B1;          % estimate weight matrix
                    C2 = (W'\B2')';     % estimate channel responses
                    
                    %C2_total(f,iter,t,i,:,:) = C2;
                    
                    % shift eegs to common center
                    n2shift = ceil(size(C2,2)/2);
                    for ii=1:size(C2,1)
                        [~, shiftInd] = min(abs(posBins-tstl(ii)));
                        C2(ii,:) = wshift('1D', C2(ii,:), shiftInd-n2shift-1);
                    end
                    
                    %tf_total(f,iter,t,i,:) = mean(C2,1); % average shifted channel responses
                    
                    % Jordan - after training on both, test on conditions
                    % seperately
                    
                    if i == 1 || i == 2
                        B2_cond1 = dt_cond1(cond1_blockNum==i,:);
                        C2_cond1 = (W'\B2_cond1')'; % test
                        %C2_total_cond1(f,iter,t,cond1Total_i,:,:) = C2_cond1;
                        
                        %shift conditions
                        n2shift = ceil(size(C2_cond1,2)/2);
                        for ii=1:size(C2_cond1,1)
                            [~, shiftInd] = min(abs(posBins-cond1_tstl(ii)));
                            C2_cond1(ii,:) = wshift('1D', C2_cond1(ii,:), shiftInd-n2shift-1);
                        end
                        
                        % origin
                        %tf_total_cond1(f,iter,t,cond1Total_i,:) = mean(C2_cond1,1);
                        tf_total_cond1(iter,trLoop,teLoop,cond1Total_i,:) = mean(C2_cond1,1);
                        
                        %C2_total_cond1Shifted(f,iter,t,cond1Total_i,:,:) = C2_cond1;
                        % Jordan added for JOCN reviewer to look at
                        % lateralized channel responses
                        
                        cond1Total_i = cond1Total_i + 1;
                        
%                         %% JORDAN
%                         tf_total_rest(iter,trLoop,teLoop,restTotal_i,:) = mean(C2_rest,1);
%                         restTotal_i = restTotal_i + 1;
%                         %% JORDAN
                        
                    elseif i == 3 || i == 4
                        
                        B2_cond2 = dt_cond2(cond2_blockNum==i,:);
                        C2_cond2 = (W'\B2_cond2')'; % test
                        %C2_total_cond2(f,iter,t,cond2Total_i,:,:) = C2_cond2;
                        
                        
                        n2shift = ceil(size(C2_cond2,2)/2);
                        for ii=1:size(C2_cond2,1)
                            [~, shiftInd] = min(abs(posBins-cond2_tstl(ii)));
                            C2_cond2(ii,:) = wshift('1D', C2_cond2(ii,:), shiftInd-n2shift-1);
                        end
                        
                        
                        %tf_total_cond2 (f,iter,t,cond2Total_i,:)  = mean(C2_cond2,1);
                        tf_total_cond2(iter,trLoop,teLoop,cond2Total_i,:) = mean(C2_cond2,1);

                        
                        %C2_total_cond2Shifted(f,iter,t,cond2Total_i,:,:) = C2_cond2;
                        
                        cond2Total_i = cond2Total_i + 1;
                        
                        
                    elseif i == 5 || i == 6
                        
                        B2_cond3 = dt_cond3(cond3_blockNum==i,:);
                        C2_cond3 = (W'\B2_cond3')'; % test
                        %C2_total_cond3(f,iter,t,cond3Total_i,:,:) = C2_cond3;
                        
                        
                        n2shift = ceil(size(C2_cond3,2)/2);
                        for ii=1:size(C2_cond3,1)
                            [~, shiftInd] = min(abs(posBins-cond3_tstl(ii)));
                            C2_cond3(ii,:) = wshift('1D', C2_cond3(ii,:), shiftInd-n2shift-1);
                        end
                        
                        
                        %tf_total_cond3 (f,iter,t,cond3Total_i,:)  = mean(C2_cond3,1);
                        tf_total_cond3(iter,trLoop,teLoop,cond3Total_i,:) = mean(C2_cond3,1);

                        
                        %C2_total_cond3Shifted(f,iter,t,cond3Total_i,:,:) = C2_cond3;
                        
                        cond3Total_i = cond3Total_i + 1;
                        
                    elseif i == 7 || i == 8
                        
                        B2_cond4 = dt_cond4(cond4_blockNum==i,:);
                        C2_cond4 = (W'\B2_cond4')'; % test
                        %C2_total_cond4(f,iter,t,cond4Total_i,:,:) = C2_cond4;
                        
                        
                        n2shift = ceil(size(C2_cond4,2)/2);
                        for ii=1:size(C2_cond4,1)
                            [~, shiftInd] = min(abs(posBins-cond4_tstl(ii)));
                            C2_cond4(ii,:) = wshift('1D', C2_cond4(ii,:), shiftInd-n2shift-1);
                        end
                        
                        
                        %tf_total_cond4 (f,iter,t,cond4Total_i,:)  = mean(C2_cond4,1);
                        tf_total_cond4(iter,trLoop,teLoop,cond4Total_i,:) = mean(C2_cond4,1);

                        
                        %C2_total_cond4Shifted(f,iter,t,cond4Total_i,:,:) = C2_cond4;
                        
                        cond4Total_i = cond4Total_i + 1;
                        
                    end
                    
                    
                    
                    
                    
                end
                % Average weights across blocks add to matrix (Mary
                % MacLean, added)
                %allWeightsTotal(:,:,:,t,f,iter) = WeightsTotal;
            end
        end
    end
    toc % stop timing the frequency loop
end

% average over iterations and blocks to reduce size of saved file
%tf_total = squeeze(mean(mean(tf_total,2),4));

tf_total_cond1 = squeeze(mean(mean(tf_total_cond1,1),4));
tf_total_cond2 =  squeeze(mean(mean(tf_total_cond2,1),4));
tf_total_cond3 = squeeze(mean(mean(tf_total_cond3,1),4));
tf_total_cond4 =  squeeze(mean(mean(tf_total_cond4,1),4));


%% compute slopes
thisX = 0:45:180; % WTF use real angular values

for iCross=1:4
    
    if      iCross==1; theseData = tf_total_cond1;
    elseif  iCross==2; theseData = tf_total_cond2;
    elseif  iCross==3; theseData = tf_total_cond3;
    elseif  iCross==4; theseData = tf_total_cond4;
    end
    
    for trSamps=1:40
        for teSamps=1:40
            %dat=rDat.total(iCross,trSamps,teSamps,:);
            dat = theseData(trSamps,teSamps,:);
            x = thisX;
            d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
            
            % for AUC?
            %fit = trapz(d);
            
            % typical slope calc.
            fit = polyfit(x,d,1);
            
            % store slopes
            em.rSl.evoked(iCross,trSamps,teSamps)=fit(1);
            
%             if      permReal==1; 
%             elseif  permReal==2; em.pSl.total(iCross,trSamps,teSamps)=fit(1);
%             end
        end
    end
end
 

% tf_total_cond1 = squeeze(mean(mean(tf_total_cond1,2),4));
% tf_total_cond2 =  squeeze(mean(mean(tf_total_cond2,2),4));
% tf_total_cond3 = squeeze(mean(mean(tf_total_cond3,2),4));
% tf_total_cond4 =  squeeze(mean(mean(tf_total_cond4,2),4));

cd([root '/Analysis_Scripts'])
% save data
if processForIEM==1
    
    fName = [root, '/IEM_Results_TT_Within_Fixed_Gen_Evoked/',sprintf('sj%02d_fixed_IEM.mat',sjNum)];
    em.C2.total = C2_total;
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
    
    em.tfs_cond1.evoked = tf_total_cond1;
    em.tfs_cond2.evoked = tf_total_cond2;
    em.tfs_cond3.evoked = tf_total_cond3;
    em.tfs_cond4.evoked = tf_total_cond4;
    
    %parsave(fName,em,minCnt,nElectrodes);
    save(fName,'em','minCnt','nElectrodes','-v7');
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


clear 
close all
end