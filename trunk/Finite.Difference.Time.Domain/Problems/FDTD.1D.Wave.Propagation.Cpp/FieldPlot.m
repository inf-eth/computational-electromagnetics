% Plots field snapshots as a movie. Taken from appendices of
% Understanding FDTD, J. B Schneider.
basename = './FieldData/Ez';
y_min = -1;
y_max = 1;
simTime = 255;
frame = 1;
figure(1);
i = 0;
data = 1:1024;
figure(1)
plothandle = plot (data);
axis ([0 length(data) y_min y_max]);
while i < simTime
    filename = sprintf ('%s%d.fdt', basename, frame);
    fid = fopen (filename, 'rt');
    if fid == -1
        return;
    end
    data = fread (fid, 1024, 'double');
    fclose (fid);
    %plot (data)
    drawnow
    set(plothandle, 'YData', data(:,1));
    %reel (frame) = getframe;
    frame = frame+1;
    i = i+1;
    pause(0.01)
end
%movie (reel, 1);
