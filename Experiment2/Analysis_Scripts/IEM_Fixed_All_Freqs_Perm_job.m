%{
Purpose:run fixed IEM script (permuted version) in parallel
Author: Tom Bullock (using Jordan's Fixed IEM script adaptation)
Date: 01.29.22

%}

%% select subjects
subjects = [1:2,4:17,20,21:24,27:28]; % kicked out subjects 03 (bad EEG) and 19 (bad BEH)

% run in serial for debugging (0) or parallel (1)
runInParallel=1;

if runInParallel
    s = parcluster();
    s.ResourceTemplate = '--ntasks-per-node=6 --mem=65536';
    job = createJob(s);
end

for iSub = 1:length(subjects)
    sjNum = subjects(iSub);
    if runInParallel
        createTask(job,@IEM_Fixed_All_Freqs_Perm,0,{sjNum});
    else
        IEM_Fixed_All_Freqs_Perm(sjNum) 
    end
end

if runInParallel
    submit(job)
    %wait(job)
end