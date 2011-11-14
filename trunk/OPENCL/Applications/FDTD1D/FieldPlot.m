% Plots field snapshots as a movie. Taken from appendices of
% Understanding FDTD, J. B Schneider.
basename = './FieldData/Ez';
y_min = -1;
y_max = 1;
simTime = 1024;
frame = 1;
figure(1);
i = 0;
while i < simTime
    filename = sprintf ('%s%d.fdt', basename, frame);
    fid = fopen (filename, 'rt');
    if fid == -1
        return;
    end
    data = fscanf (fid, '%f');
    fclose (fid);
    plot (data)
    axis ([0 length(data) y_min y_max]);
    reel (frame) = getframe;
    frame = frame+1;
    i = i+1;
end
%movie (reel, 1);
