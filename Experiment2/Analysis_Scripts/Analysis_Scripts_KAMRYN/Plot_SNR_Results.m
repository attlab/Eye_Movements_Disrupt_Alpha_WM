%{
Plot_SNR_Results
Author: Tom Bullock
Date: 05.30.20

%}

clear
close all


sourceDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';
destDir = '/home/waldrop/Desktop/WTF_EYE/Plots';


% plot baseline corrected (1) and uncorrected (2) versions
for iPlot=1:2
    
    load([sourceDir '/' 'SNR_Data.mat'])

    if iPlot==1
        % baseline correct
        snrMat = snrMat - mean(snrMat(:,:,:,77:128),4);
        these_yLims = [-.15,.07];
    else
        these_yLims = [0,.3];
    end
    
    h=figure('OuterPosition',[450   283   904   535]);
    
    % plot
    for iCond=1:4
        
        if iCond==1; thisLineStyle='-'; thisColor = 'r'; % S/F
        elseif iCond==2; thisLineStyle='-'; thisColor = 'b'; % C/F
        elseif iCond==3; thisLineStyle='-'; thisColor = 'g'; % S/M
        elseif iCond==4; thisLineStyle='-'; thisColor = 'm'; % C/M
        end
        
        %     plot(linspace(-200,2000,564),squeeze(mean(mean(snrMat(:,iCond,:,77:640),1),3)),...
        %         'color',thisColor,...
        %         'linewidth',3); hold on
        
        theseDataMean = squeeze(mean(mean(snrMat(:,iCond,:,77:640),1),3));
        theseDataSEM = squeeze(std(mean(snrMat(:,iCond,:,77:640),3),0,1))./sqrt(size(snrMat,1));
        
        shadedErrorBar(linspace(-200,2000,564),theseDataMean,theseDataSEM,...
            {'color',thisColor,...
            'linewidth',3}); hold on
        
        
    end
    
    set(gca,...
        'xlim',[-200,2000],...
        'XTick',[-200,0,500,1000,1500,2000],...
        'XTickLabel',{'-.2 ',' 0 ',' .5 ','  1 ',' 1.5','  2 '},...
        'ylim',these_yLims,...
        'linewidth',1.5,...
        'fontsize',28,...
        'box','off');
    
    pbaspect([2,1,1]);
    
    
    
    % stim offset line
    for iLine=1:5
        
        if      iLine==1; thisX=0;
        elseif  iLine==2; thisX=250;
        elseif  iLine==3; thisX=500;
        elseif  iLine==4; thisX=1000;
        elseif  iLine==5; thisX=1500;
        end
        
        line([thisX,thisX],these_yLims,...
            'linewidth',3,...
            'color','k',...
            'linestyle','--');
        
    end
    
    %legend('Spatial/Fix','Color/Fix','Spatial/Move','Color/Move','location','northwest')
    if iPlot==1
        saveas(h,[destDir '/' 'SNR_Plot_Final_BL_Corrected.eps'],'epsc')
    else
        saveas(h,[destDir '/' 'SNR_Plot_Final_BL_Uncorrected.eps'],'epsc')
    end
    
    %ylabel('SNR')
    %xlabel('Time (ms')
    %title('SNR - Total Alpha vs. All Non-Alpha Freqs')
    
end

