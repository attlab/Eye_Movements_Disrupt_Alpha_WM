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
       
       
       for iXY=1:2
           for freqs=1:30
               
               % apply butterworth filter (3rd order, bandpass)
               sampleRate = 500;
               [z1,p1] = butter(3, [freqs, freqs+1]./(sampleRate/2),'bandpass');
               
               if iXY==1
                   data = double(ex);
               else
                   data = double(ey);
               end
               
               for y = 1:size(data,1)
                   try
                       dataFilt1 = filtfilt(z1,p1,data(y,:));
                       bandEYE(y,:) = dataFilt1; % TRIED FLIPPING THIS
                   catch
                       %disp(['Trial ' num2str(y) ' contains NANS'])
                   end
               end
               
               % apply hilbert to bandpassed data
               allEye(iSub,iCond,freqs,iXY,:) = nanmean(abs(hilbert(bandEYE')).^2,2);  % mean dir here?
               
               clear bandEYE
               
           end
       end
       
   end
 
end

% save data
save([destDir '/' 'EYE_Spectral_Analysis_Hilbert.mat'],'allEye')



%%%

% load data
load([destDir '/' 'EYE_Spectral_Analysis_Hilbert.mat'])

% run t-tests against zero
for iCond=1:size(allEye,2)
    for iFreq=1:size(allEye,3)
        for iXY=1:size(allEye,4)
            for iTime=1:size(allEye,5)
                all_stats(iCond,iFreq,iXY,iTime) = ttest(allEye(:,iCond,iFreq,iXY,iTime));
            end
        end
    end
end

% generate heatmaps
h=figure('units','normalized','outerposition',[0,0,1,1]);
cnt=0;
for iXY=1:2
    for iCond=1:4
        
        if     iCond==1; thisTitle = 'Spatial Fix';
        elseif iCond==2; thisTitle = 'Color Fix';
        elseif iCond==3; thisTitle = 'Spatial Move';
        elseif iCond==4; thisTitle = 'Color Move';
        end
        
        cnt=cnt+1;
        subplot(2,4,cnt);
        imagesc(squeeze(mean(allEye(:,iCond,:,iXY,:),1)))
        set(gca,'ydir','normal','xTick',[1,114,372,620,867,1114],'xticklabel',[-250,0,500,1000,1500,2000]);
        pbaspect([1,1,1])
        title(thisTitle,'fontsize',20);
        colormap('jet');
        xlabel('Time (ms)')
        ylabel('Freq (Hz)')
        cbar
    end
end

% plot stats map
h=figure('units','normalized','outerposition',[0,0,1,1]);
cnt=0;
for iCond=1:4
    for iXY=1:2
        cnt=cnt+1;
        subplot(2,4,cnt);
        imagesc(squeeze(all_stats(iCond,:,iXY,:)),[0,1])
        cbar
    end
end


saveas(h,[plotDir '/' 'EYE_Time_Freq_Plot.jpg'],'jpeg')