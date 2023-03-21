events_to_remove = [
    94, 162, 183, 239, 270, 357, 411, 436, 501, 505, 515, 579, 602, 606, 609, 635, 700, ...
    705, 706, 711, 724, 733, 785, 812, 874, 928, 962, 969, 976, 994, 1008, 1016, 1025, ...
    1039, 1046, 1053, 1078, 1098]%, 1148, 1158, 1178, 1183, 1186];
   
allBehStructOriginal = allBehStruct;
allBehStruct(events_to_remove)=[];


EEG = pop_select(EEG, 'notrial',236);


EEG = pop_select(EEG, 'notrial',[542:623]);

EEG = pop_select(EEG, 'notrial', 568);
EEG = pop_select(EEG, 'notrial',590);
EEG = pop_select(EEG, 'notrial',593);
EEG = pop_select(EEG, 'notrial',595);
EEG = pop_select(EEG, 'notrial',597);



% check that EEG and trial data match by checking event codes
clear consistencyCheck
for i=1:length(EEG.epoch)
    consistencyCheck(i,:) = [EEG.epoch(i).eventtype{1},allBehStruct(i).locTrigger,EEG.epoch(i).eventtype{1}-allBehStruct(i).locTrigger];
end
