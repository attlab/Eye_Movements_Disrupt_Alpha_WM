%{
Beh_Merge_Blocks
Purpose: merge blocks
Author: Tom Bullock, UCSB Attention Lab
Date last modified: 11.30.21

%}
clear
close all
%function Beh_Merge_Blocks(subjects)

%which subject(s)
subjects = [27:28]; %[1:7]

% set dirs
rDir='/home/waldrop/Desktop/SCAMS';
behDir=[rDir '/BEH_Raw'];
behDirProcessed=[rDir '/' 'BEH_Data_Processed'];

% rename raw files
Beh_Rename_Files(behDir)

for iSub=1:length(subjects)
    sjNum=subjects(iSub);
    
    clear allTrialData
    
    for iSession=1:2
        
        % find CB order
        load([behDir '/' sprintf('sj%02d_se%02d_bl01_mc01_ec02.mat',sjNum,iSession)]);
        cbOrder=trialInfo(1).counterbalancingOrder;
        clear trialInfo
        
        % add exception
        if sjNum==20; cbOrder = [4,3,2,1]; end
        if sjNum==11; cbOrder = [3,4,1,2]; end
        if sjNum==27; cbOrder = [4,3,2,1]; end
        
        for iCond=cbOrder
            
            clear allTrialData
            

            
%             if      iCond==1; mc=1;ec=2;
%             elseif  iCond==2; mc=1;ec=6;
%             elseif  iCond==3; mc=2;ec=2;
%             elseif  iCond==4; mc=2;ec=6;
%             end
            
            if      iCond==1; mc=1;ec=2;
            elseif  iCond==2; mc=2;ec=2;
            elseif  iCond==3; mc=1;ec=6;
            elseif  iCond==4; mc=2;ec=6;
            end
            
            clear demographicInfo
            
            nBlocks = 4;
            
            if sjNum==20 && iCond==2; nBlocks=3; end % add exception for sj20
            
            for iBlock=1:nBlocks
                
                clear trialInfo 
                
                load([behDir '/' sprintf('sj%02d_se%02d_bl%02d_mc%02d_ec%02d.mat',sjNum,iSession,iBlock,mc,ec) ])
                
                if iBlock == 1
                    allTrialData(1:size(trialInfo,2)) = trialInfo;
                else
                    allTrialData(length(allTrialData)+1:length(allTrialData)+length(trialInfo))=trialInfo;
                end
                
            end
            
            % compile into master structure
            masterStruct(iSession,iCond).allTrialData=allTrialData;
            
        end
        
    end
    
    save([behDirProcessed '/' sprintf('sj%02d_allBeh.mat',sjNum)],'demographicInfo','cbOrder','masterStruct','-v7.3')
    
end

%return