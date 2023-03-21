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
subjects = [1:7,9:14,16:20,22:27,31];

% set dirs
rDir = '/home/waldrop/Desktop/WTF_EYE';
sourceDirReal = [rDir '/' 'IEM_Results_TT_Within_Fixed' ];
sourceDirPerm = [rDir '/' 'IEM_Results_TT_Within_Fixed_Perm' ];
%destDir = [rDir '\\' 'IEM_Slopes_TT_Within'];
destDirCompiled = [rDir '/' 'Data_Compiled'];

% subjects 
%subjects = [4,5];

% specify x values
thisX = 0:45:180; 

% sub loop
for iSub=1:length(subjects)
    sjNum = subjects(iSub);
    
    % load the "within" real data 
    real=load([sourceDirReal '/' sprintf('sj%02d_fixed_IEM.mat',sjNum)],'em');
    
    % load the "within" perm data 
    perm=load([sourceDirPerm '/' sprintf('sj%02d_fixed_IEM.mat',sjNum)],'em');

    
    
    for iCond=1:4
        
        if iCond==1
            em_within.tfs.total = real.em.tfs_cond1.total;
            em_within.tfs_perm.total = perm.em.tfs_cond1.total;
        elseif iCond==2
            em_within.tfs.total = real.em.tfs_cond2.total;
            em_within.tfs_perm.total = perm.em.tfs_cond2.total;
        elseif iCond==3
            em_within.tfs.total = real.em.tfs_cond3.total;
            em_within.tfs_perm.total = perm.em.tfs_cond3.total;
        elseif iCond==4
            em_within.tfs.total = real.em.tfs_cond4.total;
            em_within.tfs_perm.total = perm.em.tfs_cond4.total;
        end
        
        % real total data 
        f=1:2;
            for samp = 1:size(em_within.tfs.total,1)
                dat = squeeze(em_within.tfs.total(samp,:));
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                %em_within.rSl.total(f,samp)= fit(1);
                %allSlopes_real_total(iSub,iCond,f,samp) = fit(1);
                rSl_total(iSub,iCond,f,samp) = fit(1);
            end
        
        
%         % real evoked data
%         for f = 1:size(em_within.tfs.evoked,1)
%             for samp = 1:size(em_within.tfs.evoked,2)
%                 dat = squeeze(em_within.tfs.evoked(f,samp,:));
%                 x = thisX;
%                 d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
%                 fit = polyfit(x,d,1);
%                 %em_within.rSl.evoked(f,samp)= fit(1);
%                 %allSlopes_real_evoked(iSub,iCond,f,samp) = fit(1);
%                 rSl_evoked(iSub,iCond,f,samp) = fit(1);
%             end
%         end
        
        
        % perm total data
        for f = 1:2%:size(em_within.tfs_perm.total,1)
            for samp = 1:size(em_within.tfs_perm.total,1)
                dat = squeeze(em_within.tfs_perm.total(samp,:));
                x = thisX;
                d = [dat(1),mean([dat(2),dat(8)]),mean([dat(3),dat(7)]),mean([dat(4),dat(6)]),dat(5)];
                fit = polyfit(x,d,1);
                %em_within.pSl.total(f,samp)= fit(1);
                %allSlopes_perm_total(iSub,iCond,f,samp) = fit(1);
                pSl_total(iSub,iCond,f,samp) = fit(1);
            end
        end
        
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





% freqs loop (alpha, theta)
for d=[1,3]
    for f=1
        for t=1:640
            dumMem = [zeros(25,1), ones(25,1), zeros(25,1), ones(25,1)];
            dumEyes = [zeros(25,1), zeros(25,1), ones(25,1), ones(25,1)];
            
            % create column vector of subjects
            sjNums = 1:size(rSl_total,1);
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
save([destDirCompiled '/' 'IEM_Slopes_Within_Fixed.mat'],...
    'rSl_total',...
    'pSl_total',...
    'rSl_total_LF',...
    'pSl_total_LF','-v7.3')
        
  
%return
