%{
Purpose: run IEM script in parallel
Author: Tom Bullock
Date: 05.08.20
Date updated: 12.02.21

%}

% select subjects
subjects = [1:2,4:17,20,21:24,27:28]; % kicked out subjects 03 (bad EEG) and 19 (bad BEH)

% run in serial for debugging (0) or parallel (1)
runInParallel=1;

if runInParallel
    s = parcluster();
    s.ResourceTemplate = '--ntasks-per-node=6 --mem=65536';
    job = createJob(s);
end

for iSpatialOrColor=1 %[1=spatial IEM, 2=color IEM (control)]
    for iSub = 1:length(subjects)
        sjNum = subjects(iSub);
        if runInParallel
            createTask(job,@IEM_Fixed_Real_No_AAR,0,{sjNum,iSpatialOrColor});
        else
            IEM_Fixed_Real_No_AAR(sjNum,iSpatialOrColor)
        end
    end
end

if runInParallel
    submit(job)
end