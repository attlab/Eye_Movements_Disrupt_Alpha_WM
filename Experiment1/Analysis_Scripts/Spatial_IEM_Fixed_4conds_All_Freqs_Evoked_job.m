%{
Purpose: run IEM script in parallel
Author: Tom Bullock
Date: 05.08.20

%}

%% select subjects
%subjects = [1:10,13,16,19:24];
subjects = [1:7,9:14,16:20,22:27,31];



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
        createTask(job,@Spatial_IEM_Fixed_4conds_All_Freqs_Evoked,0,{sjNum});
    else
        Spatial_IEM_Fixed_4conds_All_Freqs_Evoked(sjNum) 
    end
end

if runInParallel
    submit(job)
end