clear all;

load('funwave_results/case_50m_2_nw.mat','-mat');
fresults=results;
load('matts_results/analytical_nw_case_50m_2_0.01tmp.mat','-mat');
mresults=results;



%Find how many snapshots are available
nshots=length(fresults.snapshot);

tt=[];
for k=1:length(mresults.snapshot)
    tt=[tt mresults.snapshot{k}.time];
end

f=figure(1);
% vidObj = VideoWriter('newfile_wout.avi');
% vidObj.FrameRate = 3;
% vidObj.Quality=100;
% open(vidObj);
% set(gca,'nextplot','replacechildren');

for i=10:nshots
    
    for k=1:length(mresults.snapshot)-1
        if( (mresults.snapshot{k}.time<=fresults.snapshot{i}.time)&& ...
                (fresults.snapshot{i}.time<mresults.snapshot{k+1}.time) )
            
            mx=-mresults.snapshot{k}.x;
            meta=mresults.snapshot{k}.eta;
            
            mx2=-mresults.snapshot{k+1}.x;
            meta2=mresults.snapshot{k+1}.eta;
            break
        end
    end
    
    
    %Plots the bathymetry along the bay
    p1=plot(fresults.x, fresults.b, '-g');
    hold on
    
    %Plots the FUNWAVE computed wave profile
    p2=plot(fresults.x, fresults.snapshot{i}.eta,'-b');
    
    %Plots the semi-analytic computed wave profile
    p3=plot(mx, meta,'-r');
    p3=plot(mx2, meta2,'-m');
    hold off
    
    
    %Tweak some graphics
    fresults.max_runup=max(fresults.snapshot{1}.eta);
    for lp=2:size(fresults.snapshot,2)
        fresults.max_runup=max(fresults.max_runup,max(fresults.snapshot{lp}.eta));
    end
    axis([-30 500 -0.1 fresults.max_runup*1.1])
    xlabel('Distance from the shore, meters.');
    ylabel('Wave height, meters.')
    
    set(gca,'Fontsize',12)
    
    title(['Current time = ',num2str(fresults.snapshot{i}.time),' seconds.']);
    legend([p1 p2 p3],{'Bathymetry','FUNWAVE','Semi-analytic'})
    pause(0.01)
    
    if(i==1), pause(1), end
    
%     currFrame = getframe(f);
%     writeVideo(vidObj,currFrame)
end
% 
% close(vidObj);
