%{
EYE_Import_High_Res_Into_Beh_EXP1_job
Author: Tom Bullock
Date: 10.19.22

%}

subjects = [1:7,9,10,11,12,16,17,19,20,22:27,31];

for iSub=1:length(subjects)
    sjNum=subjects(iSub);
    EYE_Import_High_Res_Into_Beh_EXP1(sjNum)
end