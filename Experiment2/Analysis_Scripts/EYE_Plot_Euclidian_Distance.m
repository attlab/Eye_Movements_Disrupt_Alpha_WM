%{
Author: Tom Bullock
Date: 10.19.22
Notes: Plot Euclidian Error for Eye-Movement trials

To Do: Convert to degrees visual angle from pixels

%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/SCAMS/Data_Compiled';
destDir = '/home/waldrop/Desktop/SCAMS/Plots';

% load data
load([sourceDir '/' 'EYE_Euclidian_Distance_Master.mat'])

for iSub = 1:length(euclid_dist) 
    for iCond=1:2
        euclid_median_all(iSub,iCond,:) = median(euclid_dist(iSub,iCond).all_trials);
    end
end


% plot
h=figure('OuterPosition',[673   733   562   267]);
for iCond=1:2 
    
    if iCond==1; thisColor = 'g';
    elseif iCond==2; thisColor = 'm';
    end
    
    thisMean = squeeze(mean(euclid_median_all(:,iCond,:),1));
    thisSEM = squeeze(std(euclid_median_all(:,iCond,:),0,1)./sqrt(size(euclid_median_all,1)));
    
    % convert all points to visual angles
    for i=1:length(thisMean)
       thisMeanVisAngle(i) = visAngleCalculate_pix(thisMean(i),120); 
       thisSEMVisAngle(i) = visAngleCalculate_pix(thisSEM(i),120);
    end
    
    
    
    %plot(-499:2000,squeeze(mean(euclid_median_all(:,iCond,:),1))); hold on
    shadedErrorBar(-499:2000,thisMeanVisAngle,thisSEMVisAngle,{'color',thisColor}); hold on
           
        
end

% add lines
for iLine = 1:4
    if iLine==1; thisTime = 0;
    elseif iLine==2; thisTime = 250;
    elseif iLine==3; thisTime = 500;
    elseif iLine==4; thisTime = 1500;
    end
    line([thisTime,thisTime],[0,4],'LineWidth',2,'Color','k','LineStyle','--');
end


set(gca,...
    'box','off',...
    'LineWidth',1.5,...
    'FontSize',16,...
    'XLim',[-500,2000],...
    'XTick',[-500,0,500,1000,1500,2000],...
    'XTickLabel',[-0.5,0,0.5,1.0,1.5,2.0],...
    'YLim',[0,4]);

pbaspect([3,1,1]);

saveas(h,[destDir '/' 'SCAMS_Euclidian_Distance_Time_Series.eps'],'epsc')
    
    