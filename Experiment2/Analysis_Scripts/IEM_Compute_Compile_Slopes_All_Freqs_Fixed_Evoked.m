%{
IEM_Compute_Slopes_Within
Purpose: calculate slope values for single freq bandpass IEM analyses

Original Author:
Joshua J. Foster
joshua.james.foster@gmail.com
University of Chicago
August 17, 2015

Modified by Tom Bullock
UCSB Attention Lab

UPDATED AGAIN 05.02.20 BY TOM
%}

%function IEM_Compute_Compile_Slopes_Within(subjects)

%if nargin==0

%end

clear
close all

subjects = [1:2,4:17,20,21:24,27:28]; % kicked out subjects 03 (bad EEG) and 19 (bad BEH)

% set dirs
rDir = '/home/waldrop/Desktop/SCAMS';
sourceDirReal = [rDir '/' 'IEM_Results_TT_Within_Fixed_All_Freqs_Evoked' ];
sourceDirPerm = [rDir '/' 'IEM_Results_TT_Within_Fixed_All_Freqs_Evoked_Perm' ];

%destDir = [rDir '\\' 'IEM_Slopes_TT_Within'];
destDirCompiled = [rDir '/' 'Data_Compiled'];

% specify x values
thisX = 0:45:180; 

% sub loop
for iSub=1:length(subjects)
    sjNum = subjects(iSub)
    
    % load the "real" data only
    allReal = load([sourceDirReal '/' sprintf('sj%02d_fixed_IEM.mat',sjNum)],'em');
    
    % load the "perm" data only
    allPerm = load([sourceDirPerm '/' sprintf('sj%02d_fixed_IEM.mat',sjNum)],'em');

    
    
    for iCond=1:4
       
        if      iCond==1; theseDataReal = allReal.em.tfs_cond1.evoked; theseDataPerm = allPerm.em.tfs_cond1.evoked;
        elseif  iCond==2; theseDataReal = allReal.em.tfs_cond2.evoked; theseDataPerm = allPerm.em.tfs_cond2.evoked;
        elseif  iCond==3; theseDataReal = allReal.em.tfs_cond3.evoked; theseDataPerm = allPerm.em.tfs_cond3.evoked;
        elseif  iCond==4; theseDataReal = allReal.em.tfs_cond4.evoked; theseDataPerm = allPerm.em.tfs_cond4.evoked;
        end
        
        
        
        % real total data
        for f = 1:size(theseDataReal,1)
            for samp = 1:size(theseDataReal,2)
                dat = squeeze(theseDataReal(f,samp,:));
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                %em_within.rSl.total(f,samp)= fit(1);
                %allSlopes_real_total(iSub,iCond,f,samp) = fit(1);
                rSl_evoked(iSub,iCond,f,samp) = fit(1); 
            end
        end
        
%         % real evoked data
%         for f = 1:size(em_within.tfs.evoked,1)
%             for samp = 1:size(em_within.tfs.evoked,2)
%                 dat = squeeze(em_within.tfs.evoked(f,samp,:));
%                 x = thisX;
%                 d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
%                 fit = polyfit(x,d,1);
%                 em_within.rSl.evoked(f,samp)= fit(1);
%                 allSlopes_real_evoked(iSub,iCond,f,samp) = fit(1);
%                 rSl_evoked(iSub,iCond,f,samp) = fit(1);
%             end
%         end
        

%         % perm total data
%         for f = 1:size(em_within.tfs_perm.total,1)
%             for samp = 1:size(em_within.tfs_perm.total,2)
%                 dat = squeeze(em_within.tfs_perm.total(f,samp,:));
%                 x = thisX;
%                 d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
%                 fit = polyfit(x,d,1);
%                 %em_within.pSl.total(f,samp)= fit(1);
%                 %allSlopes_perm_total(iSub,iCond,f,samp) = fit(1);
%                 pSl_total(iSub,iCond,f,samp) = fit(1);
%             end
%         end
%         
%         % perm total data
        for f = 1:size(theseDataPerm,1)
            for samp = 1:size(theseDataPerm,2)
                dat = squeeze(theseDataPerm(f,samp,:));
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                %em_within.rSl.total(f,samp)= fit(1);
                %allSlopes_real_total(iSub,iCond,f,samp) = fit(1);
                pSl_evoked(iSub,iCond,f,samp) = fit(1);
            end
        end
        
        
        
        
        
%         
%         % perm evoked data
%         for f = 1:size(em_within.tfs_perm.evoked,1)
%             for samp = 1:size(em_within.tfs_perm.evoked,2)
%                 dat = squeeze(em_within.tfs_perm.evoked(f,samp,:));
%                 x = thisX;
%                 d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
%                 fit = polyfit(x,d,1);
%                 %em_within.pSl.evoked(f,samp)= fit(1);
%                 %allSlopes_perm_evoked(iSub,iCond,f,samp) = fit(1);
%                 pSl_evoked(iSub,iCond,f,samp) = fit(1);
%             end
%         end
        
        clear em_within
        
    end
    
    
end




% convert into R "long" format

nSubs = length(subjects);

% freqs loop (alpha, theta)
for d=[2,4];%1:4
    for f=1:30%40;%size(em_within.tfs.total,1)
        for t=1:40
            dumMem = [zeros(nSubs,1), ones(nSubs,1), zeros(nSubs,1), ones(nSubs,1)];
            dumEyes = [zeros(nSubs,1), zeros(nSubs,1), ones(nSubs,1), ones(nSubs,1)];
            
            % create column vector of subjects
            sjNums = 1:size(rSl_evoked,1);
            sjNums = [sjNums, sjNums, sjNums, sjNums]';
            
            % isolate data
            if      d==1; theseData = rSl_total(:,:,f,t);
            elseif  d==2; theseData = rSl_evoked(:,:,f,t);
            elseif  d==3; theseData = pSl_total(:,:,f,t);
            elseif  d==4; theseData = pSl_evoked(:,:,f,t);
            end

            % write out data
            if      d==1; rSl_total_LF(:,:,f,t) = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
            elseif  d==2; rSl_evoked_LF(:,:,f,t) = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
            elseif  d==3; pSl_total_LF(:,:,f,t) = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
            elseif  d==4; pSl_evoked_LF(:,:,f,t) = [theseData(:), sjNums, dumMem(:), dumEyes(:)];
            end
            
        end
    end
end




% save all long and regular format data
save([destDirCompiled '/' 'IEM_Slopes_All_Freqs_Fixed_Evoked.mat'],...
    'rSl_evoked',...
    'pSl_evoked',...
    'rSl_evoked_LF',...
    'pSl_evoked_LF','-v7.3')
        

% do some quick plots
for iPlot=1:4
    
    subplot(1,4,iPlot)
    imagesc(squeeze(mean(rSl_evoked(:,iPlot,:,:),1)));
    pbaspect([1,2,1])
end

  
%return
