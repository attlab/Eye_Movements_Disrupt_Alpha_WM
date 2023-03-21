%{
Compile Slopes from Real and Perm Generalizations to run through R (BF
T-tests)
Author: Tom Bullock
Date: 02.10.21
%}

clear
close all

% subjects 
subjects = [1:2,4:17,20,21:24,27:28]; % kicked out subjects 03 (bad EEG) and 19 (bad BEH)
%subjects = [1,2,4:17,19:24,27,28];

% set dirs
sourceDirReal = '/home/waldrop/Desktop/SCAMS/IEM_Results_TT_Within_Fixed_Gen';
sourceDirPerm = '/home/waldrop/Desktop/SCAMS/IEM_Results_TT_Within_Fixed_Gen_Perm';
destDir = '/home/waldrop/Desktop/SCAMS/Data_Compiled';

for iSub=1:length(subjects)
    
    sjNum = subjects(iSub);
    
    % load real data
    load([sourceDirReal '/' sprintf('sj%02d_fixed_IEM.mat',sjNum)])
    
    all_rSl(iSub,:,:,:) = em.rSl.total;
    
    clear em
    
    % load perm data
    load([sourceDirPerm '/' sprintf('sj%02d_fixed_IEM.mat',sjNum)])
    
    all_pSl(iSub,:,:,:) = em.pSl.total;
    
    clear em
    
end

save([destDir '/' 'Cross_TT_Slope_MASTER_Fixed.mat'],'all_pSl','all_rSl','-v7.3')

