%{
PLOT_Quick_Look_Real_IEM_Data_Single_Subs
%}

sourceDir = '/home/waldrop/Desktop/SCAMS/IEM_Results_TT_Within_Fixed';

subject = 6;

load([sourceDir '/' sprintf('sj%02d_fixed_IEM.mat',subject)])

for iPlot=1:4
    
   subplot(1,4,iPlot)
   
   if iPlot==1
       d=em.tfs_cond1.total;
   elseif iPlot==2
       d=em.tfs_cond2.total;
   elseif iPlot==3
       d=em.tfs_cond3.total; 
   elseif iPlot==4
       d=em.tfs_cond4.total;
   end
   
   imagesc(d,[-.2,1]);
   cbar
    
end