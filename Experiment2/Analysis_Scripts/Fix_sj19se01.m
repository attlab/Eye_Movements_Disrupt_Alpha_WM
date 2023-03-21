% attempt to fix sj19
for i=1:length(EEG.epoch)
    
   % get eeg events
   eegEventList(i,1) = EEG.epoch(i).eventtype{1}; 
   
   % get beh data
   eegEventList(i,2) = allBehStruct(i).locTrigger;
   
end

allBehStruct(68) = [];
allBehStruct(71) = [];
allBehStruct(73) = [];
allBehStruct(93) = [];
allBehStruct(119) = [];
allBehStruct(120) = [];
allBehStruct(208) = [];
allBehStruct(225) = [];
allBehStruct(254) = [];
allBehStruct(256) = [];
allBehStruct(272) = [];
allBehStruct(272) = [];
allBehStruct(276) = [];
allBehStruct(278) = [];
allBehStruct(285) = [];
allBehStruct(288) = [];
allBehStruct(307) = [];
allBehStruct(316) = [];
allBehStruct(344) = [];
allBehStruct(345) = [];
allBehStruct(356) = [];
allBehStruct(369) = [];
allBehStruct(369) = [];
allBehStruct(401) = [];
allBehStruct(428) = [];
allBehStruct(429) = [];

allBehStruct(453:481) = [];
allBehStruct(457) = [];
allBehStruct(486) = [];
allBehStruct(490) = [];
allBehStruct(505) = [];
allBehStruct(522) = [];
allBehStruct(525) = [];
allBehStruct(548) = [];
allBehStruct(576) = [];
allBehStruct(581) = [];
allBehStruct(587) = [];
allBehStruct(587) = [];
allBehStruct(599) = [];
allBehStruct(602) = [];
allBehStruct(616) = [];
allBehStruct(626) = [];
allBehStruct(628) = [];
allBehStruct(641) = [];
allBehStruct(648) = [];
allBehStruct(648) = [];
allBehStruct(661) = [];
allBehStruct(663) = [];
allBehStruct(675) = [];
allBehStruct(677) = [];
allBehStruct(680) = [];
allBehStruct(686) = [];
allBehStruct(725) = [];
allBehStruct(726) = [];
allBehStruct(726) = [];
allBehStruct(731) = [];

allBehStruct(762) = [];
allBehStruct(764) = [];
allBehStruct(769) = [];
allBehStruct(773) = [];
allBehStruct(776) = [];
allBehStruct(800) = [];
allBehStruct(806) = [];
allBehStruct(818) = [];
allBehStruct(824) = [];

allBehStruct(836) = [];
allBehStruct(851) = [];
allBehStruct(854) = [];
allBehStruct(858) = [];
allBehStruct(860) = [];

allBehStruct(866) = [];
allBehStruct(866) = [];

allBehStruct(869) = [];
allBehStruct(869) = [];
allBehStruct(870) = [];
allBehStruct(894) = [];
allBehStruct(897) = [];
allBehStruct(933) = [];

allBehStruct(938) = [];
allBehStruct(948) = [];
allBehStruct(950) = [];
allBehStruct(951) = [];
allBehStruct(958) = [];
allBehStruct(959) = [];
allBehStruct(962) = [];

allBehStruct(983) = [];

allBehStruct(1010) = [];
allBehStruct(1040) = [];
allBehStruct(1055) = [];

allBehStruct(1110) = [];

