% Plots field snapshots as a movie. Taken from appendices of
% Understanding FDTD, J. B Schneider.
% Input files are assumed to be in binary format.
basename = './FieldData/Ez';
y_min = -1;
y_max = 1;
simTime = 18;
size = [1024 1];
frame = 1;
figure(1);
i = 0;
while i < simTime
    filename = sprintf ('%s%d.fdt', basename, frame);
    fid = fopen (filename);
    if fid == -1
        return;
    end
    data = fread (fid, size, 'double', 0, 'n'); 
    fclose (fid);
    plot (data)
    axis ([0 length(data) y_min y_max]);
    reel (frame) = getframe;
    frame = frame+1;
    i = i+1;
end
%movie (reel, 1);
