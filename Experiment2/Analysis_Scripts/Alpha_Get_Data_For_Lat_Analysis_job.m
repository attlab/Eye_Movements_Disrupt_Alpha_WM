%{
Purpose: run IEM script in parallel
Author: Tom Bullock
Date: 05.08.20

%}

% select subjects
subjects = [1:2,4:17,19:24,27,28];


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
        createTask(job,@Alpha_Get_Data_For_Lat_Analysis,0,{sjNum});
    else
        Alpha_Get_Data_For_Lat_Analysis(sjNum) 
    end
end

if runInParallel
    submit(job)
end