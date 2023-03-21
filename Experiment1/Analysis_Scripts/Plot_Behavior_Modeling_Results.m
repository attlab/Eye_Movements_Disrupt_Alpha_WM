%{
Plot_Beh_Data_Results
Author: Tom Bullock
Date: 12.06.19
%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';
destDir = '/home/waldrop/Desktop/WTF_EYE/Plots';

% load data
load([sourceDir '/' 'Modelling_Data.mat'])

% generate plots
h=figure('Units','normalized','OuterPosition',[0,0,1,1]);
for iPlot=1:2
    
   if       iPlot==1; theseData = modelSD; thisTitle = 'Precision'; thisYlim = [4,20]; thisYtick = [4,8,12,16,20];
   elseif   iPlot==2; theseData = modelGuess; thisTitle = 'Guess Rate'; thisYlim = [0,.2]; thisYtick = [0,.04,.08,.12,.16,.20];
   end
   subplot(1,2,iPlot)
   
   % plot bars
   for i=1:4
       if i==1; thisColor = [255,0,0];
       elseif i==2; thisColor = [0,128,255];
       elseif i==3; thisColor = [0,204,0];
       elseif i==4; thisColor = [204,0,204];
       end
       bar(i,mean(theseData(:,i),1), 'FaceColor',thisColor./255); hold on
   end
   
   % plot individual data points using plotSpread
   plotSpread(theseData,'distributionMarkers',{'.'},'distributionColors',{'k'});
   set(findall(1,'type','line','color','k'),'markerSize',24) %Change marker size
   
   
   %bar(mean(theseData,1)); hold on
   errorbar(1.25:1:4.25,mean(theseData,1),std(theseData,0,1)/sqrt(size(theseData,1)),'k.',...
       'LineWidth',2.5,...
       'CapSize',0);
   
   set(gca,'ylim',thisYlim,'ytick',thisYtick,'xlim',[.5,4.5],'linewidth',1.5)%
   
   thisMean = mean(theseData,1);
   thisMean = std(theseData,0,1)/sqrt(size(theseData,1));
   
   %set(gca,xticklabels,{'1','2','3','4'})

   xticklabelsNEW = {' '};%{'S/F','S/M','C/F','C/M'};
   
   set(gca,'xticklabels',xticklabelsNEW,'XTick',[])
   
   %title(thisTitle);
   box('off')
   set(gca,'fontsize',40)
   
   pbaspect([1,1,1]);
   
end

saveas(h,[destDir '/' 'Behavior_For_MS.eps'],'epsc')

%h=figure;
%errorbar(mean(modelSD,1),std(modelSD,1)\\sqrt(size(modelSD,1)),'o')

% save plot
%saveas(h,'ALL_SUBS_BEH_SD.fig','fig')