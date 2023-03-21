%{
Purpose: run IEM script in parallel
Author: Tom Bullock
Date: 05.08.20

%}

%% select subjects
%subjects = [1:10,13,16,19:24];
subjects =13;% [1:7,9:14,16:20,22:27,31];

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
        createTask(job,@EEG_Compute_Spectra,0,{sjNum});
    else
        EEG_Compute_Spectra(sjNum) 
    end
end

if runInParallel
    submit(job)
end