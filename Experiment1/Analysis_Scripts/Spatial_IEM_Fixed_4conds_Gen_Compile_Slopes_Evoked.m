%{
Compile Slopes from Real and Perm Generalizations to run through R (BF
T-tests)
Author: Tom Bullock
Date: 02.10.21
%}

clear
close all

% subjects 
subjects = [1:7,9:14,16:20,22:27,31];

% set dirs
sourceDirReal = '/home/waldrop/Desktop/WTF_EYE/IEM_Results_TT_Within_Fixed_Gen_Evoked';
sourceDirPerm = '/home/waldrop/Desktop/WTF_EYE/IEM_Results_TT_Within_Fixed_Gen_Perm_Evoked';
destDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';

for iSub=1:length(subjects)
    
    sjNum = subjects(iSub);
    
    % load real data
    load([sourceDirReal '/' sprintf('sj%02d_fixed_IEM.mat',sjNum)])
    
    all_rSl(iSub,:,:,:) = em.rSl.evoked;
    
    clear em
    
    % load perm data
    load([sourceDirPerm '/' sprintf('sj%02d_fixed_IEM.mat',sjNum)])
    
    all_pSl(iSub,:,:,:) = em.pSl.evoked;
    
    clear em
    
end

save([destDir '/' 'Cross_TT_Slope_MASTER_Fixed_Evoked.mat'],'all_pSl','all_rSl','-v7.3')

