%{
Eye_Frequency_FFT_Plot
Author: Tom Bullock
Date: 11.17.21
%}

for iXY = 1:2
    
    subplot(1,2,iXY);
    
    if iXY==1
        thisSpectra = allSpectra_ex;
    else
        thisSpectra = allSpectra_ey;
    end
    
    for iCond=1:4
        
        if       iCond==1; thisColor = 'r';
        elseif   iCond==2; thisColor = 'b';
        elseif   iCond==3; thisColor = 'g';
        elseif   iCond==4; thisColor = 'm';
        end
        
        thisMean = mean(thisSpectra(:,iCond,:),1);
        thisSEM = std(thisSpectra(:,iCond,:),0,1)./sqrt(size(thisSpectra,1));
        
        shadedErrorBar(1:30,thisMean,thisSEM,{'color',thisColor}); hold on
        
        %   plot(squeeze(mean(allSpectra(:,iCond,:),1)),'color',thisColor); hold on
        
    end
    
end