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
sourceDir = '/home/waldrop/Desktop/WTF_EYE/Beh_Data_Processed_Eye_Imported';
destDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';

% set subs
subjects = [1:7,9,10,11,12,16,17,19,20,22:27,31];
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
        
        % need to figure out which eye position is which
        
        for iTrial=1:length(data)
            
            % get coords for fixation dot
            xFix = data(iTrial).eyeLocX;
            yFix = data(iTrial).eyeLocY;
            
            % screen center coords
            centerX=1024/2;
            centerY=768/2;
            
            % get fixation dot new "move" pixel locations
            fixLoc1 = [(centerX + xFix), (centerY + yFix)];
            fixLoc2 = [(centerX - xFix), (centerY - yFix)];
            
            % get eye data
            ex = data(iTrial).eyeHighRes.gx(2,:);
            ey = data(iTrial).eyeHighRes.gy(2,:);
            
            % get eye movements for each movement window
            times = -500:2000;
            theseTimes = [500,1000];
            theseTimesIdx = (find(times==theseTimes(1)):find(times==theseTimes(2)));
            
            ex1 = ex(theseTimesIdx);
            ey1 = ey(theseTimesIdx);
            
            theseTimes = [1000,1500];
            theseTimesIdx = (find(times==theseTimes(1)):find(times==theseTimes(2)));
            
            ex2 = ex(theseTimesIdx);
            ey2 = ey(theseTimesIdx);
            
            % get means over em1 and em2 windows
            mean_em1 = [mean(ex1),mean(ey1)];
            mean_em2 = [mean(ex2),mean(ey2)];
            
            % get euclidian distances of em1/em2 from fixLoc1/fixLoc2
            for j=1:length(ex1)
                
                % compare em1 to both pairs of fixLocs
                e1(j) = sqrt((fixLoc1(1)-ex1(j))^2 + (fixLoc1(2)-ey1(j))^2);
                e2(j) = sqrt((fixLoc2(1)-ex1(j))^2 + (fixLoc2(2)-ey1(j))^2);
                
                % compare em2 to both paris of fixLocs
                e3(j) = sqrt((fixLoc1(1)-ex2(j))^2 + (fixLoc1(2)-ey2(j))^2);
                e4(j) = sqrt((fixLoc2(1)-ex2(j))^2 + (fixLoc2(2)-ey2(j))^2);
                
            end
            
            % find min euclidian distance of EMs to fixLoc1 and fixLoc2
            [min_em1_fixLoc1, min_em1_fixLoc1_idx] = min(e1);
            [min_em1_fixLoc2, min_em1_fixLoc2_idx] = min(e2);
            [min_em2_fixLoc1, min_em2_fixLoc1_idx] = min(e3);
            [min_em2_fixLoc2, min_em2_fixLoc2_idx] = min(e4);
            
            %
            if min_em1_fixLoc1<min_em1_fixLoc2
                firstFixLoc = 1;
                fixLoc1_color = 'g';
                this_em1_ed = min_em1_fixLoc1;
                this_em1_pix = [ex1(min_em1_fixLoc1_idx),ey1(min_em1_fixLoc1_idx)];
                this_em1_ed_time = min_em1_fixLoc1_idx;
            else
                firstFixLoc = 2;
                fixLoc1_color = 'b';
                this_em1_ed = min_em1_fixLoc2;
                this_em1_pix = [ex1(min_em1_fixLoc2_idx),ey1(min_em1_fixLoc2_idx)];
                this_em1_ed_time = min_em1_fixLoc2_idx;
            end
            
            if firstFixLoc == 1
                fixLoc2_color = 'b';
                this_em2_ed = min_em2_fixLoc2;
                this_em2_pix = [ex2(min_em2_fixLoc2_idx),ey2(min_em2_fixLoc2_idx)];
                this_em2_ed_time = min_em2_fixLoc2_idx;
            else
                fixLoc2_color = 'g';
                this_em2_ed = min_em2_fixLoc1;
                this_em2_pix = [ex2(min_em2_fixLoc1_idx),ey2(min_em2_fixLoc1_idx)];
                this_em2_ed_time = min_em2_fixLoc1_idx;
            end
            
            % save euclidian distance errors to a struct
            euclidErrorStruct(iCond).em1(iTrial) = this_em1_ed;
            euclidErrorStruct(iCond).em2(iTrial) = this_em2_ed;
            
            %euclidErrorStruct(iCond).em1_idx(iTrial) = this_em1_ed_time;
            %euclidErrorStruct(iCond).em2_idx(iTrial) = this_em2_ed_time;
            
            
            for eyeMovement = 1:2
                
                if      firstFixLoc==1 && eyeMovement==1
                    cuePosXY =  fixLoc1;
                elseif  firstFixLoc==2 && eyeMovement==1
                    cuePosXY = fixLoc2;
                elseif firstFixLoc==2 && eyeMovement==2
                    cuePosXY = fixLoc1;
                elseif firstFixLoc==1 && eyeMovement==2
                    cuePosXY = fixLoc2;
                end
                   
                
                
            end
            
            % GIVING UP FOR NOW COZ THIS PROBS WON'T WORK
            
            
            
            
            
            
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