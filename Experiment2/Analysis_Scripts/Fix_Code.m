% attempt to fix sj11 se01
for i=1:length(EEG.epoch)
   % get eeg events
   eegEventList(i,1) = EEG.epoch(i).eventtype{1}; 
   % get beh data
   eegEventList(i,2) = allBehStruct(i).locTrigger;
end