function [FigWidth,FigHeight] = Cfd2lcsPlotRoms(fname,varname,cmin,cmax,FigWidth,FigHeight)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot a contour plot of a cfd2lcs variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
%fname:  Filename
%varname:  Variable name
%cmin: Min contour value
%cmax: Max contour value
%OUTPUTS:
%FigWidth:  Width of the figure
%%FigHeight:  Height of the figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

%Load the data:
x = h5read(fname,'/CFD_GRID/-X');
y = h5read(fname,'/CFD_GRID/-Y');
z = h5read(fname,'/CFD_GRID/-Z');
flag = h5read(fname,'/FLAG');
ftle = h5read(fname,varname);

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
meanX=mean(x');
meanY=mean(y);

%Compute the correct scale of the figure the first time:
if(FigHeight<=0)
    FigHeight = round(FigWidth*(max(max(y))-min(min(y))) / (max(max(x))-min(min(x)))) ;
end

%Zero the ftle where flag is zero:
ftle(flag>0)=0.0;

%The Figure:
ss = get(0,'screensize'); %The screen size
sw = ss(3);
sh = ss(4);
fh = figure('Position',[ sw-FigWidth sh-FigHeight FigWidth FigHeight]);
levels = linspace(cmin,cmax,10);
contourf(meanX,meanY,ftle',levels,'LineStyle','none','LevelStepMode','manual');
shading flat;
colormap hot;
colormap(flipud(colormap));
caxis([cmin cmax]);
hold on;
[c,h]=contour(meanX,meanY,flag','k','LineWidth',2); %Boundary
hold on;
contour(meanX,meanY,deglat',[0:1:180],'LineStyle',':'); %Latitude grid
contour(meanX,meanY,deglon',[0:1:180],'LineStyle',':'); %Longitude grid
set(gca,'xtick',[],'ytick',[])
end