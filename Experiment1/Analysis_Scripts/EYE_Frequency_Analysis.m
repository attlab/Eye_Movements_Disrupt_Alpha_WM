%{
EYE_Frequency_Analyses
Author: Tom Bullock
Date: 11.17.21

%}

clear
close all

% set dirs
sourceDir = '/home/waldrop/Desktop/WTF_EYE/EYE/Eye_Sync_Final';
destDir = '/home/waldrop/Desktop/WTF_EYE/Data_Compiled';
plotDir = '/home/waldrop/Desktop/WTF_EYE/Plots';

% set subs
subjects = [1:3,5:7,9:12,16,17,19,20,22,24:27,31];% all look normal (others either missing data or look crazy)

for iSub = 1:length(subjects)
    
    iSub
    
    %figure;
    
   sjNum = subjects(iSub);
   
   % load data
   load([sourceDir '/' sprintf('sj%02d_eye_beh_sync_pro.mat',sjNum)])
   
   % combine two sessions into single mat
   for iCond=1:4
       
       ex = [];
       ey = [];
       spectra_ex = [];
       spectra_ey = [];
       ex_trial = [];
       ey_trial = [];
       tmp = [];
       
       % combine sessions
       ex = [eyeProcessed(1,iCond).ex;eyeProcessed(2,iCond).ex];
       ey = [eyeProcessed(1,iCond).ey;eyeProcessed(2,iCond).ey];
       
       % run fft [do separate x/y for now, then try euclid dist from fix]      
       NFFT = 500;
       L = 1114; 
       Fs = 500;
       
       for iTrial=1:size(ex,1)
           ex_trial = ex(iTrial,:);
           ey_trial = ey(iTrial,:);
           spectra_ex(iTrial,:) = fft(ex_trial,NFFT)/L;
           spectra_ey(iTrial,:) = fft(ey_trial,NFFT)/L;
       end
       
       if       iCond==1; thisColor = 'r';
       elseif   iCond==2; thisColor = 'b';
       elseif   iCond==3; thisColor = 'g';
       elseif   iCond==4; thisColor = 'm';
       end
       
%          % induced ?
%         tmp = plot(nanmean(abs(spectra(:,2:31)).^2),'color',thisColor); hold on        
%         title(num2str(sjNum))
        
        allSpectra_ex(iSub,iCond,:) = nanmean(abs(spectra_ex(:,2:31)).^2);
        allSpectra_ey(iSub,iCond,:) = nanmean(abs(spectra_ey(:,2:31)).^2);
        
       
        % evoked ?
       %tmp = plot(abs(nanmean(spectra(:,2:31)).^2),'color',thisColor); hold on
       
       %set(gca,'YScale','log')
       
       % get freqs
       myFreqs = Fs/2*linspace(0,1,NFFT/2+1); %frequencies for the single sided power spectrum
       
   end
 
end

save([destDir '/' 'EYE_Spectral_Analysis.mat'],'myFreqs','allSpectra_ex','allSpectra_ey')

% plot
h=figure;%('units','normalized','outerposition',[0,0,1,1]);
for iXY = 1:2
    
    subplot(1,2,iXY);
    
    if iXY==1
        thisSpectra = allSpectra_ex;
        thisTitle = 'Eye X';
    else
        thisSpectra = allSpectra_ey;
        thisTitle = 'Eye Y';
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
    
    
    xlabel('Freq (Hz)')
    ylabel('Power (uV^2)')
    title(thisTitle,'Fontsize',20)
    
end

saveas(h,[plotDir '/' 'EYE_Frequency_Plot.jpeg'],'jpeg')