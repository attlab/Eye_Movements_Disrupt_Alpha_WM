function EYE_Import_High_Res_Into_Beh(sjNum)
%{
EYE_Import_High_Res_Into_Beh
Author: Tom Bullock (UCSB Attention Lab)
Date: 10.18.22


%}

%clear
%close all

 % set subject
%sjNum=7;

% set dirs
sourceDirBeh = '/home/waldrop/Desktop/SCAMS/BEH_Data_Processed';
sourceDirEye = '/home/waldrop/Desktop/SCAMS/EYE_Raw';
destDir = '/home/waldrop/Desktop/SCAMS/Beh_Data_Processed_Eye_Imported';

% load beh data
load([sourceDirBeh '/' sprintf('sj%02d_allBeh.mat',sjNum)])
disp(['Loading Beh File: ' sprintf('sj%02d_allBeh.mat',sjNum)])

% loop through all the sessions and factorial combinations and grab
% eyetracker data using initial target timestamps
upsampledDataPresent = 0;
for iSession=1:2
    for iMemCond=1:2
        for iEyeCond=[2,6]
            
            if iMemCond==1 && iEyeCond==2
                condCode=1;
            elseif iMemCond==2 && iEyeCond==2
                condCode=2;
            elseif iMemCond==1 && iEyeCond==6
                condCode=3;
            elseif iMemCond==2 && iEyeCond==6
                condCode=4;
            end
            
            % load eye data
            load([sourceDirEye '/' sprintf('s%d_%02d%d%d.mat',sjNum,iSession,iMemCond,iEyeCond)])
            disp(['Loading Eye File ' sprintf('s%d_%02d%d%d.mat',sjNum,iSession,iMemCond,iEyeCond)])  
            
            % upsample eye data if accidentally recorded at 250 Hz instead of 1000 Hz
            if edfStruct.FSAMPLE.time(2) - edfStruct.FSAMPLE.time(1) == 4
                disp('Upsample Data')
                edfStruct.FSAMPLE.time = repelem(edfStruct.FSAMPLE.time,1,4);
                edfStruct.FSAMPLE.pa = repelem(edfStruct.FSAMPLE.pa,1,4);
                edfStruct.FSAMPLE.gx = repelem(edfStruct.FSAMPLE.gx,1,4);
                edfStruct.FSAMPLE.gy = repelem(edfStruct.FSAMPLE.gy,1,4);
                upsampledDataPresent = 1;
            end
            
            % grab high res eye data
            for iTrial=1:length(masterStruct(iSession,condCode).allTrialData)
                if masterStruct(iSession,condCode).allTrialData(iTrial).brokenTrial~=1
                    startTime = masterStruct(iSession,condCode).allTrialData(iTrial).targEye.elSampleTime(1);
                    
                    if ~ismember(sjNum,[7,8,9,12]) % upsampled subs have to be treated differently
                        trialTimeIdx = find(edfStruct.FSAMPLE.time==startTime-500) : find(edfStruct.FSAMPLE.time==(startTime+1999));
                    else
                        [val,startIdx] = min(abs(double(edfStruct.FSAMPLE.time) - (startTime-500)));
                        [val,endIdx] = min(abs(double(edfStruct.FSAMPLE.time) - (startTime+1999)));
                        trialTimeIdx = startIdx:endIdx;
                        if length(trialTimeIdx>2500); trialTimeIdx = trialTimeIdx(1:2500);end
                    end
                    
                    masterStruct(iSession,condCode).allTrialData(iTrial).eyeHighRes.time = edfStruct.FSAMPLE.time(trialTimeIdx);
                    masterStruct(iSession,condCode).allTrialData(iTrial).eyeHighRes.pa = edfStruct.FSAMPLE.pa(:,trialTimeIdx);
                    masterStruct(iSession,condCode).allTrialData(iTrial).eyeHighRes.gx = edfStruct.FSAMPLE.gx(:,trialTimeIdx);
                    masterStruct(iSession,condCode).allTrialData(iTrial).eyeHighRes.gy = edfStruct.FSAMPLE.gy(:,trialTimeIdx);
                else
                    masterStruct(iSession,condCode).allTrialData(iTrial).eyeHighRes = [];
                end
            end
            
        end
    end
end

% save data to new folder
save([destDir '/' sprintf('sj%02d_allBeh.mat',sjNum)],'masterStruct','upsampledDataPresent','-v7.3')

return