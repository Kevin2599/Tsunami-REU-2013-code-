t=[88,189,205,220,300];
freq=10;
freq2=100;
p1=plot(fresults.x, fresults.b, '-g');
hold on
for i=1:5
    p2=plot(fresults.x, fresults.snapshot{t(i)}.eta,'-b');
    p3=plot(-mresults.snapshot{t(i)}.x, mresults.snapshot{t(i)}.eta,'-r');
    if i==2
        p4=plot(fresults.x(1:freq2:end), fresults.snapshot{t(i)}.eta(1:freq2:end),'^b');
        p5=plot(-mresults.snapshot{t(i)}.x(1:freq:end), mresults.snapshot{t(i)}.eta(1:freq:end),'^r');
    end
    if i==3
        plot(fresults.x(1:freq2:end), fresults.snapshot{t(i)}.eta(1:freq2:end),'ob');
        plot(-mresults.snapshot{t(i)}.x(1:freq:end), mresults.snapshot{t(i)}.eta(1:freq:end),'or');
    end
    if i==4
        p6=plot(fresults.x(1:freq2:end), fresults.snapshot{t(i)}.eta(1:freq2:end),'*b');
        p7=plot(-mresults.snapshot{t(i)}.x(1:freq:end), mresults.snapshot{t(i)}.eta(1:freq:end),'*r');
    end
    if i==5
        plot(fresults.x(1:freq2:end), fresults.snapshot{t(i)}.eta(1:freq2:end),'.b','Markersize',20);
        plot(-mresults.snapshot{t(i)}.x(1:freq:end), mresults.snapshot{t(i)}.eta(1:freq:end),'.r','Markersize',20);
    end
    
axis([-30 300 -0.1 fresults.max_runup*1.1])
xlabel('Distance from the shore, meters.');
ylabel('Wave height, meters.')
pause(2)
end
legend([p1 p2 p3 p4 p5 p6 p7],{'Bathymetry','FUNWAVE','Semi-analytic','FUNWAVE At min rundown','Semi-analytic At min rundown','FUNWAVE At max runup','Semi-analytic At max runup'})

hold off