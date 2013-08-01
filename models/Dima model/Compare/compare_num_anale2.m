clear all;

load('funwave_results/funwave.mat','-mat');
fresults=results;
load('matts_results/analytical.mat','-mat');
mresults=results;

fresults.max_runup=0.95;

%Find how many snapshots are available
nshots=length(fresults.snapshot);

tt=[];
for k=1:length(mresults.snapshot)
    tt=[tt mresults.snapshot{k}.time];
end

max_runup=-1;
for i=1:length(fresults.snapshot)
    iii=abs(fresults.snapshot{i}.eta-fresults.b)<1e-2;
    fresults.snapshot{i}.eta(iii)=NaN;
    max_runup=max(max_runup, max(fresults.snapshot{i}.eta));
end





t=[138,189,205,221,300];
t=[60 90 100 108 120 150];
freq=10;
freq2=100;
f=figure(1); clf
set(f,'Position',[100 100 800 1000])



index=1:length(t);
for i=index
    s(i)=subplot(length(index),1,i)
    p1=plot(fresults.x, fresults.b, '-g');
    hold on
    plot([0 1500], [0 0], '--','Color',[1 1 1]*0.5)
    
    plot(fresults.x, fresults.snapshot{t(i)}.eta,'-b');
    plot(-mresults.snapshot{t(i)}.x, mresults.snapshot{t(i)}.eta,'-r');
    
    plot(fresults.x(1:freq2:end), fresults.snapshot{t(i)}.eta(1:freq2:end),'ob');
    plot(-mresults.snapshot{t(i)}.x(1:freq:end), mresults.snapshot{t(i)}.eta(1:freq:end),'vr');
    
%     plot([-30 600], [max_runup max_runup],'-r');
    
    p4=plot([-100 -100], [-1 1],'-ob');
    p5=plot([-100 -100], [-1 1],'-vr');
    hold off
    f=1;
    if(max(fresults.snapshot{t(i)}.eta)<max_runup/2), f=2; end
        
    axis([-25 400 -max_runup*1.1/f max_runup*1.1/f])
    text(310,0.78*max_runup/f, ['Time = ',num2str(t(i)),' sec'], 'FontSize', 13)
    
    if(i==length(index)), xlabel('Distance from the shore, m', 'FontSize', 13),
    else set(s(i),'XTickLabels',[]); end;
    
    if(i==1) legend([p4 p5], {'FUNWAVE','Semi-analytic'},...
                             'Location','NorthOutside',...
                             'Orientation','Horizontal'); end
    
    set(gca, 'FontSize',12)
    ylabel('Wave height, m', 'FontSize', 13)
    pause(.02)

end

for i=index
    set(s(i),'Position', [0.13 0.97-i/length(index)/1.1 0.775 1/length(index)/1.3])
end





return

f=figure(1);
% vidObj = VideoWriter('newfile_wout.avi');
% vidObj.FrameRate = 3;
% vidObj.Quality=100;
% open(vidObj);
% set(gca,'nextplot','replacechildren');

for i=221:1:221 %nshots
    
    for k=1:length(mresults.snapshot)-1
        if( abs(mresults.snapshot{k}.time-fresults.snapshot{i}.time)<1e-5 )
            mx=-mresults.snapshot{k}.x;
            meta=mresults.snapshot{k}.eta;
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
    hold off
    
    
    %Tweak some graphics
    axis([-30 500 -0.1 fresults.max_runup*1.1])
    xlabel('Distance from the shore, meters.');
    ylabel('Wave height, meters.')
    
    set(gca,'Fontsize',12)
    
    title(['Current time = ',num2str(fresults.snapshot{i}.time),' seconds.']);
    legend([p1 p2 p3],{'Bathymetry','FUNWAVE','Semi-analytic'})
    pause(0.1)
    
    if(i==1), pause(1), end
    
    %     currFrame = getframe(f);
    %     writeVideo(vidObj,currFrame)
end
%
% close(vidObj);
