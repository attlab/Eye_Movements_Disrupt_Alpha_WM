function SpatialIEM_TrainBoth_Generalize_Permute (subs)
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

Modified by Jordan Garrett
%}
%==========================================================================

%clear all
close all

%subs = [1:8,10];


em.nPerms =100; % how many permutations? %%%TOM EDIT THIS WAS SET TO 10!!%%%
nPerms = em.nPerms;
%% process for IEM (1=IEM)
processForIEM=1;

% setup directories
root = '/home/garrett/WTF_Bike/';
dRoot = [root, 'Data/TrainBoth/'];
eRoot = [root,'EEG/'];
hRoot = [root,'HILBERT/'];

name = '_SpatialTF_ALPHA_Gen.mat'; % name of files to be saved for IEM
freqBandpass = [8 12];

EEG = [];
eegInfo = [];

% Grab TF data file
fName = [dRoot, sprintf('sj%02d_TrainBoth_changeDect_acc',subs), name];
thisFile = load(fName);
em = thisFile.em;
thisFile = [];

em.time = -.5*1000:1000/em.Fs:1.9961*1000; %   -500:4:2000; % time points of interest

% get analysis settings from TF data file.
nChans = em.nChans;
nBins = em.nBins;
nIter = em.nIter;
%nIter=10;
nBlocks = em.nBlocks;
freqs = freqBandpass;
times = em.time;
nFreqs = size(em.frequencies,1);
%%nElectrodes = em.nElectrodes;
nSamps = length(em.time);
Fs = em.Fs;
basisSet = em.basisSet;
posBin = em.posBin;
nTrialsPerBlock = em.nTrialsPerBlock;

genSamps = em.genSamps; %generalization samples

fprintf('Subject:\t%d\n',subs)

% Grab data------------------------------------------------------------

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
trial_indx = [repmat(1,1,length(rest_eeg.newTrialInfo)), repmat(2,1,length(ex_eeg.newTrialInfo))];
trial_indx = trial_indx(~artInd);
trial_indx = trial_indx';

em.nTrials = length(posBin); 
nTrials = em.nTrials; % # of good trials


%----------------------------------------------------------------------
% Preallocate arrays
permInd = nan(nFreqs,nIter,nPerms,nBlocks,nTrialsPerBlock);
permedBins = nan(1,nTrialsPerBlock);
tf_total_rest=[];
tf_total_ex = [];

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
        blocks = em.blocks(:,iter); % grab blocks assignment for current iteration
        
        % Loop through permutations
        for perm = 1:nPerms
            tic
            
            %-----------------------------------------------------------------------------
            % Permute trial assignment within each block
            %-----------------------------------------------------------------------------
            permedPosBin = nan(size(posBin)); % preallocate permuted position bins vector
            
            %Jordan added
            perm_Trials = nan(size(trial_indx));
            pInd_Trial = [];
            
            for b = 1:nBlocks % for each block..
                pInd = randperm(nTrialsPerBlock); % create a permutation index;
                permedBins(pInd) = posBin(blocks == b); % grab block b data and permute according data according to the index
                permedPosBin(blocks == b) = permedBins; % put permuted data into permedPosBin
                permInd(f,iter,perm,b,:) = pInd; % save the permutation (permInd is saved at end of the script)
                
                %Jordan added
                pInd_Trial(pInd) = trial_indx(blocks == b);
                perm_Trials(blocks == b) = pInd_Trial;
                
            end
            
            
            %-------------------------------------------------------------------------
            
            % Average data for each position bin across blocks
            posBins = 1:nBins;
            
            blockDat_total = nan(nBins*nBlocks,nElectrodes,nSamps);  % averaged total data
            
            blockDat_total_Rest = nan(nBins*nBlocks/2,nElectrodes,nSamps);
            blockDat_total_Low = blockDat_total_Rest;
            
            labels = nan(nBins*nBlocks,1);                           % bin labels for averaged data
            blockNum = nan(nBins*nBlocks,1);                         % block numbers for averaged data
            c = nan(nBins*nBlocks,nChans);                           % predicted channel responses for averaged data
            bCnt = 1;
            
            % Jordan added
            rest_labels = [];%nan(nBins*nBlocks/2,1);
            rest_blockNum = [];%nan(size(rest_labels));
            rest_c =[];%nan(nBins*nBlocks/2,nChans);
            
            ex_labels = nan(size(rest_labels));
            ex_blockNum = nan(size(ex_labels));
            ex_c = nan(size(rest_c));
            
            rest_bCnt = 1;
            ex_bCnt = 1;
            for ii = 1:nBins
                for iii = 1:nBlocks
                    
                    blockDat_total(bCnt,:,:) = squeeze(mean(fdata_total(permedPosBin==posBins(ii) & blocks==iii,:,tois),1));
                    
                    
                    % Seperate condition data
                    if (iii == 1 || iii == 2) % Rest condition
                        
                        blockDat_total_Rest(rest_bCnt,:,:) = abs(squeeze(mean(fdata_total(permedPosBin==posBins(ii) & blocks==iii,:,tois),1)));
                        
                        
                        % Con labels for testing
                        rest_labels(rest_bCnt) = ii;
                        rest_blockNum(rest_bCnt) = iii;
                        rest_c(rest_bCnt,:) = basisSet(ii,:);
                        
                        rest_bCnt = rest_bCnt + 1;
                    elseif (iii == 3 || iii == 4) % Exercise condition
                        
                        blockDat_total_Low(ex_bCnt,:,:) = abs(squeeze(mean(fdata_total(permedPosBin==posBins(ii) & blocks==iii,:,tois),1)));
                        
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
            for trLoop=1:genSamps %divide 625 samples into 25 samples of 25 (100ms per sample)
                
                this_trainSamps = (genSamps*trLoop)-(genSamps-1):(genSamps*trLoop);% divide 640 samples into 40 samples of 16 (62.5ms per sample)
                
                trainData = squeeze(mean(blockDat_total(:,:,this_trainSamps),3));
                
                
                % testing loop
                for teLoop=1:genSamps
                    
                    this_testSamps = (genSamps*teLoop)-(genSamps-1):(genSamps*teLoop);
                    
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
                            
                            tf_total_rest(perm,iter,trLoop,teLoop,restTotal_i,:) = mean(C2_rest,1);
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
                            
                            tf_total_ex(perm,iter,trLoop,teLoop,exTotal_i,:) = mean(C2_ex,1);
                            exTotal_i = exTotal_i + 1;
                        end
                    end
                end
            end
        end
    end
end

toc % stop timing the frequency loop

% For non-parametric statistics, average over iterations and then blocks
tf_total_rest = squeeze(mean(mean(tf_total_rest,2),5));
tf_total_ex = squeeze(mean(mean(tf_total_ex,2),5));


% save data
fName = [dRoot, sprintf('sj%02d_TrainBoth_changeDect_acc_Permute',subs), name];
em.permInd = permInd;
em.permtfs.total.rest = tf_total_rest;
em.permtfs.total.ex = tf_total_ex;
em.permtfs.nBlocks = nBlocks;
em.permtfs.c=c;
save(fName,'em','nElectrodes','-v7.3');




close all
clear all

