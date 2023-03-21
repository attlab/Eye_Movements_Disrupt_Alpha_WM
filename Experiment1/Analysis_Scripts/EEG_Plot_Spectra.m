%{
EEG_Spectra_Plot
Author: Tom Bullock
Date: 11.18.21

set up to plot a comparison of post AAR vs no AAR

%}

clear
close all

for iPlot=1:2
    
    subplot(1,2,iPlot);
    
    if iPlot==1
        sourceDir = '/home/waldrop/Desktop/WTF_EYE/Spectra_Post_AAR';
    else
        sourceDir = '/home/waldrop/Desktop/WTF_EYE/Spectra_No_AAR';
    end
    destDir = '/home/waldrop/Desktop/WTF_EYE/Plots';
    
    subjects = [1:7,9:14,16:20,22:27,31];
    
    % compile spectra
    for iSub=1:length(subjects)
        sjNum = subjects(iSub);
        load([sourceDir '/' sprintf('sj%02d_spectra.mat',sjNum)])
        spectraCombined(iSub,:,:,:) = allSpectra;
    end
    
    % plot
    for iCond=1:4
        
        if       iCond==1; thisColor = 'r';
        elseif   iCond==2; thisColor = 'b';
        elseif   iCond==3; thisColor = 'g';
        elseif   iCond==4; thisColor = 'm';
        end
        
        thisMean = squeeze(mean(mean(spectraCombined(:,iCond,:,:),1),3));
        thisSEM = squeeze(mean(std(spectraCombined(:,iCond,:,:),0,1),3)./sqrt(size(spectraCombined,1)));
        
        shadedErrorBar(freqs,thisMean,thisSEM,{'color',thisColor},1); hold on
        set(gca,'ylim',[0,4]);
        
        if iPlot==1
            title('WITH AAR')
        else
            title('NO AAR')
        end
        
    end
    
end