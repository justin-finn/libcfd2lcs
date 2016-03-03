clear vars; close all;
%%%%%Create an Animation:
FigWidth = 400;
Prefix = './../res0/cfd2lcs_output/bkwdFTLE_';
Suffix = '.h5';
VarName = '/BKWD-FTLE';
MovieFile = 'FTLE.avi';
MovieTitle = 'Backward FTLE (T = 7 Days)';
FPS = 5;
StartFrame = 160;
EndFrame = 210;
ContourMin = 2e-6;
ContourMax = 4.5e-6;
Margin = 5;
%%%%%

%Open movie file:
v = VideoWriter(MovieFile);
v.FrameRate = FPS;
open(v);

%Loop through the data:
FigHeight=0;
for Frame = StartFrame:EndFrame
    %Load the contour var:
    %Zero the ContourVar where flag is zero:
    FileName = sprintf('%s%d%s',Prefix,Frame,Suffix)
    flag = h5read(FileName,'/FLAG');
    ContourVar = h5read(FileName,VarName);
    ContourVar(flag>0)=0.0;
   
    %Open Figure and set x,y,z on the first frame:
    if(Frame==StartFrame)
        %Load the Coordinates:
        x = h5read(FileName,'/GRID/GRID-X');
        y = h5read(FileName,'/GRID/GRID-Y');
        z = h5read(FileName,'/GRID/GRID-Z');
    
        %Convert to spherical, then to mercator projection
        [lon,lat,R] = cart2sph(x,y,z);
        deglon = lon*180/pi;
        deglat = lat*180/pi;
        x = lon;
        y = log(abs(tan(lat)+sec(lat)));
        scaleFactor = sec(y);
        meanSF = mean(mean(scaleFactor)); % use the mean scalefactor (a single number)
        meanR = mean(mean(R)); %use the mean scalefactor (a single number)
        x=x*meanR/meanSF; %X in Units of M
        y=y*meanR/meanSF; %Y in Units of M
        %Plot Lat/Lon on axes::
        Xcoord = mean(x').*(deglon(:,1)'./mean(x')); 
        Ycoord = mean(y).*(deglat(1,:)./mean(y)); 
               
        %Determine the axes rectangle for export.
        FigHeight = round(FigWidth*(max(max(y))-min(min(y))) / (max(max(x))-min(min(x)))) ;
        ss = get(0,'screensize'); %The screen size
        fh = figure('Position',[ ss(3)-FigWidth ss(4)-FigHeight FigWidth FigHeight]);
        contour(Xcoord,Ycoord,flag','k','LineWidth',2); %Boundary
        ax = gca;
        ax.Units = 'pixels';
        ax.XTick = [0:1:180];
        ax.YTick = [0:1:180];
        ax.XLabel.String = 'Degrees Longitude';
        ax.YLabel.String = 'Degrees Latitude';
        ax.Title.String = MovieTitle;
        drawnow
        pos = ax.Position;
        ti = ax.TightInset;
        rect = [-ti(1)-Margin, -ti(2)-Margin, pos(3)+ti(1)+ti(3)+2*Margin, pos(4)+ti(2)+ti(4)+2*Margin];
    end      

    %The Color contour:
    levels = linspace(ContourMin,ContourMax,10);
    [C,h] = contourf(Xcoord,Ycoord,ContourVar',levels,'LineStyle','none','LevelStepMode','manual');      
    shading interp;
    colormap hot;
    colormap(flipud(colormap));
    caxis([ContourMin ContourMax]);
    %Boundary contour
    hold on;
    contour(Xcoord,Ycoord,flag','k','LineWidth',2); %Boundary
    %Lat/Lon grid:
    contour(Xcoord,Ycoord,deglat',[0:1:180],'LineStyle',':');
    contour(Xcoord,Ycoord,deglon',[0:1:180],'LineStyle',':');
    
    %Set the Axis labels:        
    ax.XTick = [0:1:180];
    ax.YTick = [0:1:180];
    ax.XLabel.String = 'Degrees Longitude';
    ax.YLabel.String = 'Degrees Latitude';
    ax.Title.String = MovieTitle;
    
    hold off;
    
    %Grab this frame:
    drawnow
    frame = getframe(ax,rect); 
    writeVideo(v,frame);
end

%Close movie file
close(v);
