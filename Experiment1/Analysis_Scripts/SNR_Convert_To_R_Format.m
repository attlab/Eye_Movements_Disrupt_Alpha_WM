%{
SNR_Convert_To_R_Format
Author: Tom Bullock
Date: 06.18.21
%}

clear
close all

sourceDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';
destDir = sourceDir;

% load data
load([sourceDir '/' 'SNR_Data.mat'])

% select channels for averaging
elects = [25,26,27,29,30,62:64];

snrMat = squeeze(mean(snrMat(:,:,elects,:),3));

% loop through times
for iTime=1:640
    
    thisData = snrMat(:,:,iTime);
    
    % convert into R "long" format
    dumMem = [zeros(25,1), ones(25,1), zeros(25,1), ones(25,1)];
    dumEyes = [zeros(25,1), zeros(25,1), ones(25,1), ones(25,1)];
    
    % create column vector of subjects
    sjNums = 1:size(thisData,1);
    sjNums = [sjNums, sjNums, sjNums, sjNums]';
    
    % use colon operator function to get into long format
    model_SNR_LF(:,:,iTime) = [thisData(:), sjNums, dumMem(:), dumEyes(:)];
    
end

% save data
save([destDir '/' 'SNR_Data_For_R.mat'],'model_SNR_LF','-v7.3')