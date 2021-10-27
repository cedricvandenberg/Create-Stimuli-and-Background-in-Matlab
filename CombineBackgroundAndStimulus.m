function CombineBackgroundAndStimulus;
%Read background & stimulus pattern, scale, place stimulus at random on background
close all;
clear all;

[files,nf]=GetFileList('','STIM*.mat');
[idx,tf]=listdlg('ListString',files,'SelectionMode','single',...
    'PromptString','CLICK on a pattern name then ok','ListSize',[160 150],...
    'Name','Pattern Select');
stf=files(idx); stf=char(stf);
load(stf); 
name=fn; clear files nf idx tf stf k fn;
%PATTERN VARIABLES *****************************************************************
% name       name of design
% np         number of patch outlines
% edges{np}  patch outline x,y coordinates
% codes(np)  patch colour codes
% clrs(nc,3) RGB values for each color
%***********************************************************************************
%rescale stimulus to have a total length of about 0.1 (1/10 of background)
for p=1:np
    edg=edges{p}; x=edg(:,1); y=edg(:,2); 
    x=x/8000; y=y/8000; %total length about 800 
    edg(:,1)=x; edg(:,2)=y;
    edges{p}=edg;
end;
fprintf(1,'Using stimulus %s\n',name);
figure(1); set(gcf,'Position',[90 32 560 240]);
for p=1:np
  edg=edges{p}; x=edg(:,1); y=edg(:,2); cd=codes(p);
  patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
end; 
axis equal; axis off;
title(['Stimulus: ' name]);

[files,nf]=GetFileList('','BACK*.mat');
[idx,tf]=listdlg('ListString',files,'SelectionMode','single',...
    'PromptString','CLICK on a background file name then ok',...
    'ListSize',[200 100],'Name','Background Select');
stf=files(idx); stf=char(stf); 
load(stf); backfile=fn; clear files nf idx stf tf tt fn prompt;
%PATTERN VARIABLES *****************************************************************
% name       name of design
% np         number of patch outlines
% edges{np}  patch outline x,y coordinates
% codes(np)  patch colour codes
% clrs(nc,3) RGB values for each color
%BACKGROUND VARIABLES
% backfile   background file name
% ttl        details of background generation
% nt         number of patches in background
% xx,yy(nt)  x,y coordinates of each background patch
% patchcolor(nt,3) RGB values of each background patch
%***********************************************************************************

%get a random location between 0.2 and 0.8 on both axes
rs=0.2+0.6*rand(2,1); sx=rs(1); sy=rs(2); 

% get an orientation angle (degrees) and size multiplier
tt='Orientation & Size';
prompt={'Enter R (random) or orientation angle (deg, 0=horizontal)',...
        'Enter size multiplier (1=no change from 1/10)'};
definput={'R','1.0'};
os=inputdlg(prompt,tt,[1,46],definput); os=char(os); 
str=os(1);
mlt=str2num(os(6))/10+0.04;
k=strfind(str,'R'); if isempty(k) k=0; end; 
if k==0 
  rot=str2num(os(1));
else
  rot=360*rand;  
end;
clear k tt prompt definput os;

%find centroid of stimulus by collecting all points
cx=[]; cy=[];
for p=1:np
  edg=edges{p}; x=edg(:,1); y=edg(:,2);
  cx=[cx; x]; cy=[cy; y];
end;
[pcx,pcy]=polycenter(cx,cy); clear cx cy p edg x y; % hold on; plot(pcx,pcy,'ok');
%center on 0,0 and rotate outline to rot angle
for p=1:np
  edg=edges{p}; x=edg(:,1)-pcx; y=edg(:,2)-pcy; %move centroid to origin
  [rx,ry]=RotateXYpoints(x,y,0,0,rot); %rotate by an angle rot counterclockwise
  rx=mlt*rx; ry=mlt*ry; %scale to new size relative to about 0.1 length
  edg(:,1)=rx; edg(:,2)=ry;
  edges{p}=edg;
end;
clear rx ry x y rot;

fprintf(1,'Making background image display, takes minutes to show\n');
figure(2); set(gcf,'Position',[90 32 1226 1098]); set(gcf, 'InvertHardCopy', 'off');
patch([0 1 1 0],[0 0 1 1],[0.6 0.6 0.6],'EdgeColor','none');
hold on;
% show patches
for k=1:nt
   x=xx{k}; y=yy{k}; c=patchcolor(k,:); 
   patch(x,y,c,'EdgeColor','none');    
end;
axis equal; axis off;
title(ttl);
drawnow;

%add the stimulus to the background at sx,sy
hold on;
for p=1:np
  edg=edges{p}; x=edg(:,1)+sx; y=edg(:,2)+sy; cd=codes(p);
  patch(x,y,clrs(cd,:),'EdgeColor',clrs(cd,:)); 
end; 
axis equal; axis off; drawnow;
title(['Stimulus: ' name]);

%save the panel as a .png file
prompt={'Enter file name to save'}; tt='FILE NAME';
fn=inputdlg(prompt,tt,[1 30]); fn=char(fn); fn=[fn '.png'];
print(figure(2),fn,'-dpng','-r300');
fprintf(1,'Background saved in %s\n',fn);
close all;
end %function

function [files,nf]=GetFileList(dirn,stype)
% [files,nf]=GetFileList(dirn,stype);
% INPUT (strings in single-quotes)
%   dirn is a directory name, such as 'C:\Active\Heinsohn\'  % Must end '\'
%    to read in the current directory, use '' (empty)
%   stype is the kind of file, such as '.ttt' or '.Master.transmission'
% OUTPUT  
%   files is a list of the files (without the directory in front)
%   nf    is the number of files (can be 0)
% Note: will recognize .Master.tra with stype='.tra'
drp=[dirn '*' stype];
drl=dir(drp); [nf,f]=size(drl);
files=[];
for f=1:nf
  temp=drl(f).name;
  files=[files; cellstr(temp)];
end;
end %get file list

function [rx,ry]=RotateXYpoints(x,y,ox,oy,rot)
% [rx,ry]=RotateXYpoints(x,y,ox,oy,rot);
% Rotate a set of points a given angle around a given point
% INPUT
%  x,y    points to rotate
%  ox,oy  rotate x,y around this point
%  rot    angle of rotation (degrees)  positive is clockwise (up to 360)
rot=deg2rad(rot);             %degrees to radians
cx=x-ox; cy=y-oy;             %center on the rotating center (move to zero)
[angs,dsts]=cart2pol(cx,cy);  %angles & distances of the points
angs=angs-rot;                %rotate clockwise by angle rot (radians)
[rx,ry]=pol2cart(angs,dsts);  %convert back to cartesian
rx=rx+ox; ry=ry+oy;           %move back to original position
end %rotate xy points
