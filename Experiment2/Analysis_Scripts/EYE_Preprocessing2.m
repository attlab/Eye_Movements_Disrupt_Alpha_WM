%{
Process_Eye_Data2
Author: Tom Bullock
Date:04.30.20

Sync eye with trial data
%}

clear
close all

eyeDir = '/home/waldrop/Desktop/SCAMS/EYE_Processed1';
behDir = '/home/waldrop/Desktop/SCAMS/BEH_Data_Processed';
saveDir = '/home/waldrop/Desktop/SCAMS/Eye_Synchronized';

%%subjects = [21:24,27,28]; %[1:6,10:11,13:17]; % some of these are not working - need to look into why that is
subjects = [1:2,4:17,20,21:24,27:28]; % kicked out subjects 03 (bad EEG) and 19 (bad BEH)

% NOTE SUBJECT 20 IS MISSING A FILE - SKIP IN ANALYSIS


for iSub=1:length(subjects)
    
    try
        
        % select subject number (just one)
        sjNum=subjects(iSub)
        
        load(sprintf([behDir '/sj%02d_allBeh.mat'],sjNum))
        
        for e=1:8
            
            if      e==1; eyeFile=sprintf('s%d_0112.mat',sjNum); a=1; b=1;
            elseif  e==2; eyeFile=sprintf('s%d_0122.mat',sjNum); a=1; b=2;
            elseif  e==3; eyeFile=sprintf('s%d_0116.mat',sjNum); a=1; b=3;
            elseif  e==4; eyeFile=sprintf('s%d_0126.mat',sjNum); a=1; b=4;
            elseif  e==5; eyeFile=sprintf('s%d_0212.mat',sjNum); a=2; b=1;
            elseif  e==6; eyeFile=sprintf('s%d_0222.mat',sjNum); a=2; b=2;
            elseif  e==7; eyeFile=sprintf('s%d_0216.mat',sjNum); a=2; b=3;
            elseif  e==8; eyeFile=sprintf('s%d_0226.mat',sjNum); a=2; b=4;
            end
            
            load([eyeDir '/' eyeFile])
            
            keepGoing=1;
            i=0;
            cnt=0;
            while keepGoing
                
                i=i+1;
                
                beh=masterStruct(a,b).allTrialData(i).locTrigger;
                eye=eyeEpoch.targCode(i);
                
                if beh~=eye && eye==0
                    %disp('let it be')
                elseif beh==eye
                    %disp('correct')
                else
                    masterStruct(a,b).allTrialData(i)=[];
                    cnt=cnt+1;
                    badBehVec(cnt)=i;
                end
                
                checkEvents(i,1:2) = [masterStruct(a,b).allTrialData(i).locTrigger, eyeEpoch.targCode(i)];
                checkEvents(i,3) = checkEvents(i,1)-checkEvents(i,2);
                
                if checkEvents(i,3)~=0 || checkEvents(i,3)<200
                else
                    disp('BAD WRONG YOU MESSED UP BIG TIME!!!')
                end
                
                
                % quit while loop if we hit final value in eyeEpoch arrays
                if i==length(eyeEpoch.targCode)
                    keepGoing=0;
                end
                
            end
            
            masterStruct(a,b).eyeData = eyeEpoch;
            masterStruct(a,b).checkEvents = checkEvents;
            
            clear beh eye badBehVec checkEvents
            
        end
        
        save([saveDir '/' sprintf('sj%02d_eye_beh_sync.mat',sjNum)],'masterStruct','cbOrder','-v7.3')
        
        clear masterStruct eyeEpoch demographicInfo
        
    catch
        disp(['Error processing sj' num2str(sjNum)])
    end
    
end
