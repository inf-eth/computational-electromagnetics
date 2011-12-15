% Plots field snapshots as a movie. Taken from appendix of
% Understanding FDTD, J. B Schneider.
% Input files are assumed to be in binary format.
basename = './FieldData/Ez';

fidp = fopen ('./FieldData/Parameters.smp', 'r', 'l');
if fidp == -1
    return;
end
datap = fread (fidp, 7, 'uint');
fclose (fidp);

I = datap(1)
J = datap(2)
trez = datap(3)
xrez = datap(4)
yrez = datap(5)
simTime = datap(6)/trez-1    % Number of frames to be read. Last saved field number 
                     % If last field saved is Ez1023.fdt, maximum simTime should be 1023.
PMLw = datap(7)
J = J+2*PMLw;
                    
size = [I J];    % Spatial size or width w.
frame = 1;
reel = 0;
i = 0;
while i < simTime
    filename = sprintf ('%s%d.fdt', basename, frame);
    fid = fopen (filename, 'r', 'l');
    if fid == -1
        return;
    end
    data = fread (fid, size, 'double'); 
    figure(1);
    surf (data(1:xrez:I,1:yrez:J))
    view (0, 90)
    zlim ( [-2 2] )
    caxis([-1 1])
    
    figure(2);
    mesh (data(1:xrez:I,1:yrez:J))
    view (4, 4)
    zlim ( [-2 2] )
    caxis([-1 1])
    
%     reel (frame) = getframe;
    frame = frame+1;
    i = i+1;
    fclose (fid);
end
% movie (reel, 1);
