%{
SCAM Plot Slopes for Within IEM
Author: Tom Bullock
Date: 04.22.20

Updated to plot BF thresholded bars to indicate sig
%}

clear 
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/SCAMS/Data_Compiled';
destDir = '/home/waldrop/Desktop/SCAMS/Plots';
sourceDirBayesResults = '/home/waldrop/Desktop/WTF_EYE/R_Results';

% load data
load([sourceDir '/' 'IEM_Slopes_Within_Fixed.mat']);

% load BF results
load([sourceDirBayesResults '/' 'BF_Results_IEM_WIthin_Alpha_RvP_Fixed.mat'])

% set freq band (1=alpha,2=theta)
f=1;

%% set up plot
% line widths
observedLineWidth = 3;
nullLineWidth = 2;
axisLineWidth = 1.5;

% font sizes
axisFontSize = 18;
titleFontSize = 22;
xLabelFontSize = 24;
yLabelFontSize = 24;

% colors
critTcolor = [.5 .5 .5];
obsTcolor = [0.04 0.52 .78];

% axis limits and labels
thisSR = 256; % sample rate
xNewTick = [1, thisSR*.25, thisSR*.5, thisSR*.75, thisSR*1, thisSR*1.25, thisSR*1.5, thisSR*1.75, thisSR*2, thisSR*2.25, thisSR*2.5];
xNewTickLabel = ['-0.5';'-.25';'  0 ';' .25';' 0.5';' .75';'  1 ';'1.25';' 1.5';'1.75';'  2 '];
thisYlim = [-.0017 .008];%'auto';

% plot total slope
thisYpos = -.0006; % start position for stats plots
h=figure('Units','Normalized','OuterPosition',[0.0052    0.2026    0.4696    0.7688]); % change these numbers to change plot size
subplot(2,1,1);
for iCond=1:4
    
    if iCond==1; thisLineStyle='-'; thisColor = 'r'; % S/F
    elseif iCond==2; thisLineStyle='-'; thisColor = 'b'; % C/F
    elseif iCond==3; thisLineStyle='-'; thisColor = 'g'; % S/M
    elseif iCond==4; thisLineStyle='-'; thisColor = 'm'; % C/M
    end
    
    % plot real data
    %realLineIdx(iCond)=plot(squeeze(mean(rSl_total(:,iCond,f,:),1)),'lineStyle',thisLineStyle,'LineWidth',observedLineWidth,'color',thisColor); hold on
    
    % plot real data with error bars
    thisMean = squeeze(mean(rSl_total(:,iCond,f,:),1));
    thisSEM = squeeze(std(rSl_total(:,iCond,f,:),0,1)./sqrt(size(rSl_total,1)));
    shadedErrorBar(1:640,thisMean,thisSEM,{'lineStyle',thisLineStyle,'LineWidth',observedLineWidth,'color',thisColor},1); hold on

    if iCond==1; thisLineStyle='--'; thisColor = 'r'; % S/F
    elseif iCond==2; thisLineStyle='--'; thisColor = 'b'; % C/F
    elseif iCond==3; thisLineStyle='--'; thisColor = 'g'; % S/M
    elseif iCond==4; thisLineStyle='--'; thisColor = 'm'; % C/M
    end
    
    % plot permuted data
    %plot(squeeze(mean(pSl_total(:,iCond,f,:),1)),'lineStyle',thisLineStyle,'LineWidth',observedLineWidth,'color',thisColor); hold on
    
    % plot permuted data with error bars
    thisMean = squeeze(mean(pSl_total(:,iCond,f,:),1));
    thisSEM = squeeze(std(pSl_total(:,iCond,f,:),0,1)./sqrt(size(pSl_total,1)));
    shadedErrorBar(1:640,thisMean,thisSEM,{'lineStyle',thisLineStyle,'LineWidth',observedLineWidth,'color',thisColor},1)

    % plot stats
    bfTreshold = 3;
    thisYpos = thisYpos-.0002;
    for b=1:size(bfGrandOutput,1)
        if bfGrandOutput(b,iCond)>bfTreshold
            line([b,b+1],[thisYpos,thisYpos],...
                'LineWidth',3,...
                'LineStyle',thisLineStyle,...
                'Color',thisColor)
        end
    end
  
    
    
    % plot settings
    box('off')
    set(gca,'xTick',xNewTick,'xticklabel',xNewTickLabel,'xLim',[1 640],'yLim',thisYlim,'fontsize',axisFontSize,'LineWidth',axisLineWidth);
    vline(thisSR*.5,'k--')
    vline(thisSR*.75,'k--')
    vline(thisSR*1,'k--')
    vline(thisSR*1.5,'k--')
    vline(thisSR*2,'k--')
    xlabel('Time (s)','FontSize',xLabelFontSize)
    ylabel('Slope','FontSize',yLabelFontSize)
    box('off')
    pbaspect([3 1 1])
    %title('CRF Slopes','FontSize',titleFontSize)
end

% set legend
%h_legend = legend(realLineIdx,'Spatial/Fix','Color/Fix','Spatial/Move','Color/Move');
%set(h_legend,'FontSize',12);

% % save figure
% saveas(h,[destDir '/' 'CTF_Slope_Alpha_Within_Fixed.jpeg'],'jpeg')
% saveas(h,[destDir '/' 'CTF_Slope_Alpha_Within_Fixed.eps'],'epsc')


%% plot Bayes Factors 
%h=figure('Units','Normalized','OuterPosition',[0    0.3432    0.4832    0.4068]); % change these numbers to change plot sizea=gcf
subplot(2,1,2);
for iCond=1:4
    
    if iCond==1; thisLineStyle='-'; thisColor = 'r'; % S/F
    elseif iCond==2; thisLineStyle='-'; thisColor = 'b'; % C/F
    elseif iCond==3; thisLineStyle='-'; thisColor = 'g'; % S/M
    elseif iCond==4; thisLineStyle='-'; thisColor = 'm'; % C/M
    end
    
    % plot BFs
    
    plot(bfGrandOutput(:,iCond),'lineStyle',thisLineStyle,'LineWidth',observedLineWidth,'color',thisColor); hold on
    
    % plot settings
    box('off')
    set(gca,'xTick',xNewTick,'xticklabel',xNewTickLabel,'xLim',[1 640],'ylim',[0.2108 1.0000e+12],'yScale','log','fontsize',axisFontSize,'LineWidth',axisLineWidth);
    vline(thisSR*.5,'k--')
    vline(thisSR*.75,'k--')
    vline(thisSR*1,'k--')
    vline(thisSR*1.5,'k--')
    vline(thisSR*2,'k--')
    
    %hline(0,'k--')
    hline(1,'k--')
    hline(3,'k--')
    hline(10,'k--')
    %hline(100,'k--')
    
    xlabel('Time (s)','FontSize',xLabelFontSize)
    ylabel('Bayes Factor','FontSize',yLabelFontSize)
    box('off')
    pbaspect([3 1 1])
    
    
    %realLineIdx(iCond)=plot(squeeze(mean(rSl_total(:,iCond,f,:),1)),'lineStyle',thisLineStyle,'LineWidth',observedLineWidth,'color',thisColor); hold on
    
end

% set legend
%h_legend = legend(realLineIdx,'Spatial/Fix','Color/Fix','Spatial/Move','Color/Move');
%set(h_legend,'FontSize',12);

% save figure
saveas(h,[destDir '/' 'CTF_Slope_Alpha_Within_Fixed_Plus_BFs_Error_Bars.jpeg'],'jpeg')
saveas(h,[destDir '/' 'CTF_Slope_Alpha_Within_Fixed_Plus_BFs_Error_Bars.eps'],'epsc')