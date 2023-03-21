%{
Process_Eye_Data
Author: Tom Bullock
Date: 03.16.22

Notes: epoch around Trigger 102 and create matrix of all trials 

%}

clear
close all

sourceDir = '/home/waldrop/Desktop/SCAMS/EYE_Raw';
destDir = '/home/waldrop/Desktop/SCAMS/EYE_Processed1';

cd(sourceDir)
d=dir('s*.mat');

for f=1:length(d)
    
    try
    
    %disp(num2str(f));
    
    % load data file
    load([d(f).folder '/' d(f).name])
    
    % upsample if accidentally recorded at 250 Hz instead of 1000 Hz
    if edfStruct.FSAMPLE.time(2) - edfStruct.FSAMPLE.time(1) == 4
    
        disp('Upsample Data')
        
        edfStruct.FSAMPLE.time = repelem(edfStruct.FSAMPLE.time,1,4);
        edfStruct.FSAMPLE.pa = repelem(edfStruct.FSAMPLE.pa,1,4);
        edfStruct.FSAMPLE.gx = repelem(edfStruct.FSAMPLE.gx,1,4);
        edfStruct.FSAMPLE.gy = repelem(edfStruct.FSAMPLE.gy,1,4);
    
    
    end
    
    
    
    
    % extract triggers
    for t=1:length(edfStruct.FEVENT)  
        if length(edfStruct.FEVENT(t).message)==11 %~isempty(edfStruct.FEVENT(t).message)      
            if strcmp(edfStruct.FEVENT(t).message(1:7),'TRIGGER')       
                %if length(edfStruct.FEVENT(t).message)>10   
                    triggerInfo(t,:) = str2double(edfStruct.FEVENT(t).message(end-2:end));   
                %end   
            end   
        end      
    end
    
    % extract trigger times
    triggerTimes = double([edfStruct.FEVENT.sttime]'); % will need to chop last two "end events" triggers for consistency

    % extract sample times
    sampleTimes = double(edfStruct.FSAMPLE.time);
    xPos = double(edfStruct.FSAMPLE.gx(2,:));
    yPos = double(edfStruct.FSAMPLE.gy(2,:));
    pa = double(edfStruct.FSAMPLE.pa(2,:));
    
   
%     % create variables (original import method)
%     triggerInfo = eyeMat.event.Messages.info;
%     triggerTimes = eyeMat.event.Messages.time;
%     sampleTimes = double(eyeMat.RawEdf.FSAMPLE.time);
%     xPos = double(eyeMat.RawEdf.FSAMPLE.gx(2,:));
%     yPos = double(eyeMat.RawEdf.FSAMPLE.gy(2,:));
%     pa = double(eyeMat.RawEdf.FSAMPLE.pa(2,:));
    
    
    cnt=0;
    for i=1:length(triggerInfo)
        
        if triggerInfo(i) == 102 %strcmp(triggerInfo(i),'TRIGGER 102')
            cnt=cnt+1;
            
            % get timestamp index (102) for each trial
            [~,startTimeIdx] = min(abs(triggerTimes(i)-sampleTimes));
            
            % get timestamps, x/y positions and pupil diameter
            eyeEpoch.allTimes(cnt,:) = sampleTimes(startTimeIdx:startTimeIdx+1730); % 1750 = 3500 ms post 102 (should be enough to capture full trial?)
            eyeEpoch.allXpos(cnt,:) = xPos(startTimeIdx:startTimeIdx+1730);
            eyeEpoch.allYpos(cnt,:) = yPos(startTimeIdx:startTimeIdx+1730);
            eyeEpoch.allPA(cnt,:) = pa(startTimeIdx:startTimeIdx+1730);
            
            % get code and timestamp for target (if target present)
            %if ismember(str2double(triggerInfo{i+1}(end-2:end)),201:208)   
            if ismember(triggerInfo(i+1),201:208)
                eyeEpoch.targCode(cnt,:) = triggerInfo(i+1); %str2double(triggerInfo{i+1}(end-2:end));
                eyeEpoch.targTime(cnt,:) = triggerTimes(i+1);
            else
                eyeEpoch.targCode(cnt,:) = 0;
                eyeEpoch.targTime(cnt,:) = 0;
            end
        end
        
    end
    
    save([destDir '/' d(f).name],'eyeEpoch')
    
    clear cnt eyeEpoch eyeMat pa sampleTimes triggerInfo triggerTimes xPos yPos
    
    catch
        
        disp(['Unable to process ' d(f).name])
        
    end
        
end