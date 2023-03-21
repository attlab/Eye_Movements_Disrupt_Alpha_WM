%{
EYE_Distance_Analysis
Author: Tom Bullock
Date:10.18.22

Note Screen Resolution is 1920 (w) x 1080 (h)

Try 1500-2000 for Euclid Dist Mean

%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/SCAMS/Beh_Data_Processed_Eye_Imported';
destDir = '/home/waldrop/Desktop/SCAMS/EYE_Euclidian_New_Analysis';

% set subs
subjects = [1:2,4:17,21:24,27:28]; %sj07,08,09,12 are problem coz upsample??
%subjects = [10,11,13:17,21:24,27,28];
%subjects = 7;

% define center fixation point
centerXY = [1920/2, 1080/2];

% loop through subs
for iSub=1:length(subjects)
    sjNum = subjects(iSub)
    
    % load behavioral data with high res eye data imported
    load([sourceDir '/' sprintf('sj%02d_allBeh.mat',sjNum)])
    
    % merge sessions and remove broken trials
    %%euclid_dist = [];
    for iCond=3:4
        
        if iCond==3
            thisCond=1;
        else
            thisCond=2;
        end
        
        % merge sessions
        data = [masterStruct(1,iCond).allTrialData, masterStruct(2,iCond).allTrialData];
        
        % remove broken trials
        data(find([data.brokenTrial]==1)) = [];
        
        
        for iTrial=1:length(data)
            
            % absolute EM cue position on screen in pixels
            cuePosXY = [data(iTrial).eyeLocX, data(iTrial).eyeLocY] + centerXY;
            
            % now find euclidian distance between eye and cue?
            euclid_dist(iSub,thisCond).all_trials(iTrial,:) = sqrt(((data(iTrial).eyeHighRes.gx(2,:) - cuePosXY(1)).^2) + ((data(iTrial).eyeHighRes.gy(2,:) - cuePosXY(2)).^2));
            
            %euclid_dist(thisCond).all_trials(iTrial,:)
            
%             plot(euclid_dist); hold on
%             pause(1)
            
            clear cuePosXY
            
        end
        
        clear data
        
    end
    
    %plot(median(euclid_dist(1).all_trials(:,:))); hold on
    %plot(median(euclid_dist(2).all_trials(:,:)));
    
end

save([destDir '/' 'EYE_Euclidian_Distance_Master.mat'],'euclid_dist','subjects')


% % central fixation location
% for iTrial=1:length(masterStruct(1,3).allTrialData)
%     subplot(1,2,1);
%     if masterStruct(1,3).allTrialData(iTrial).brokenTrial~=1
%         plot(-499:2000,masterStruct(1,3).allTrialData(iTrial).eyeHighRes.gx(2,:)); hold on
%         
%     end
%     subplot(1,2,2);
%     if masterStruct(1,3).allTrialData(iTrial).brokenTrial~=1
%         plot(-499:2000,masterStruct(1,3).allTrialData(iTrial).eyeHighRes.gy(2,:)); hold on
%         pause(1)
%     end
% end