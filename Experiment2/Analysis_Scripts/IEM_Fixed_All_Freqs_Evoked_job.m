%{
Purpose: run IEM script in parallel
Author: Tom Bullock
Date: 05.08.20

%}

%% select subjects
%subjects = 27; %[1,2,4:17,19:24];
%subjects = [1:2,4:17,20,21:24,27:28]; % kicked out subjects 03 (bad EEG) and 19 (bad BEH)
subjects = 6;

% run in serial for debugging (0) or parallel (1)
runInParallel=0;

if runInParallel
    s = parcluster();
    s.ResourceTemplate = '--ntasks-per-node=6 --mem=65536';
    job = createJob(s);
end

for iSub = 1:length(subjects)
    sjNum = subjects(iSub);
    if runInParallel
        createTask(job,@IEM_Fixed_All_Freqs_Evoked,0,{sjNum});
    else
        IEM_Fixed_All_Freqs_Evoked(sjNum) 
    end
end

if runInParallel
    submit(job)
end