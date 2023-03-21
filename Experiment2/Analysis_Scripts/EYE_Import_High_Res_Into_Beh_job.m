
%subjects = [1:2,4:17,20,21:24,27:28]; % kicked out subjects 03 (bad EEG)
%and 19 (bad BEH) ALSO 20 no good

for sjNum=subjects
    EYE_Import_High_Res_Into_Beh(sjNum)
end