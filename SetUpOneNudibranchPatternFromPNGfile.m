%extract nudibranch pattern outline from PNG file and save coded loops in .mat file
clear all;
close all;
% ddir='./Outlines/';  flt=[ddir '*.png'];
[fname,path]=uigetfile('./Outlines/*.png');
fn=[path fname]; img=imread(fn);
bw=im2bw(img); bw=imcomplement(bw); bw=bwskel(bw); bw=bwmorph(bw,'diag');
cc=bwconncomp(bw);  stats=regionprops(cc,'PixelList'); np=length(stats);
patches=cell(np,1);
for p=1:np
  xy=stats(p).PixelList; ty=xy(1,1); tx=xy(1,2); P=[tx ty];
  bd=bwtraceboundary(bw,P,'N'); nN=length(bd);
  if nN<10 bd=bwtraceboundary(bw,P,'E'); end;
  x=bd(:,2); y=bd(:,1); patches{p}=[x y];  %x=columns
  % plot(x,y,'b');
end;
  figure; set(gcf,'Position',[35 60 995 565]); 
  title(fname); hold on;
  for p=1:np
    xy=patches{p}; x=xy(:,1); y=xy(:,2); plot(x,y,'k');
  end;
  axis equal; axis off; title(fname);
clear P bd bw cc ddir f fn nN p path stats tx ty x y xy;
fprintf(1,'%s outline read\n',fname);
%VARIABLES *************************************************************************
% file        name of nudibranch stylized outline
% np          number of patches (outlines)
% patches{np} patch outlines, each saved as xy(:,2)
%***********************************************************************************
nclr=inputdlg('Number of colours','Colour Classes'); nclr=str2num(char(nclr));
pmpt='Body edge/ring (1)';
for c=2:nclr pmpt=[pmpt ', colour ' num2str(c) '(' num2str(c) ')']; end;


% figure; set(gcf,'Position',[835 560 995 565]); hold on;
clf; hold on;
pts=patches; np=length(pts); zns=zeros(np,1); zns(1)=1;
for p=1:np
  xy=pts{p}; 
  plot(xy(:,1),xy(:,2),'b'); axis equal;
  x=inputdlg(pmpt,'IDENTIFY EACH ZONE EDGE'); x=str2num(char(x));
  zns(p)=x;
  switch x
    case 1 
      plot(xy(:,1),xy(:,2),'b','linewidth',2); axis equal;
      text(xy(1,1),xy(1,2),num2str(x));
    case 2
      plot(xy(:,1),xy(:,2),'k','linewidth',2);
      text(xy(1,1),xy(1,2),num2str(x));
    case 3
      plot(xy(:,1),xy(:,2),'r','linewidth',2);
      text(xy(1,1),xy(1,2),num2str(x));
    case 4
      plot(xy(:,1),xy(:,2),'g','linewidth',2);
      text(xy(1,1),xy(1,2),num2str(x));
    otherwise
      plot(xy(:,1),xy(:,2),'m','linewidth',2);   
  end;
end;
zones=zns;
clear p x xy zns; 
axis equal;


%rearrange spots in order of codes for zone map saving order
cds=zones; nb=length(cds); idx=(1:nb)';
mx=[cds idx]; mx=sortrows(mx,1); idx=mx(:,2);
pts=pts(idx); cds=mx(:,1);
outlines=pts; zcodes=cds;

ttlc{1}='black=layer 1';
ttlc{2}='blue=layer 2';
ttlc{3}='red=layer 3';
ttlc{4}='green=layer 4';
ttlc{5}='magenta=layer 5';
nc=max(cds);
ttl=ttlc{1}; for c=2:nc ttl=[ttl ', ' ttlc{c}]; end;

clf
hold on; 
for b=1:nb
  xy=outlines{b}; x=xy(:,1); y=xy(:,2); cd=zcodes(b);
  switch cd
    case 1
      plot(x,y,'k'); axis equal;
      text(x(1),y(1),num2str(cd));
    case 2
      plot(x,y,'b');
      text(x(1),y(1),num2str(cd));
    case 3
      plot(x,y,'r');
      text(x(1),y(1),num2str(cd));
    case 4
      plot(x,y,'g');
      text(x(1),y(1),num2str(cd));
    otherwise
      plot(x,y,'m');
      text(x(1),y(1),num2str(cd));
  end;
end;
axis equal;
title(ttl);
drawnow;
%get mat file name from file name
k=strfind(fname,'.'); if isempty(k) k=0; end; k=k-1;
sname=fname(1:k);
clear b c* f* idx img k nb ncds nclr mx patches pmpt pts ttl ttlc x xy y zones;  
%VARIABLES **********************************************************************
% sname        pattern name
% nc           number of colours ( = maximum of zcodes)
% np           number of patches
% outlines{np} patch outlines to make zone map
% zcodes(np)   zone codes as order of laying out the zone map (smallest size later)
%******************************************************************************** 
CONTENTS={...
    'sname        pattern name'; ...
    'nc           number of colours ( = maximum of zcodes)'; ...
    'np           number of patches'; ...
    'outlines{np} patch outlines to make zone maps'; ...
    'zcodes(np)   zone codes as order of laying out the zone map (smaller later)'};
save(sname);
fprintf(1,'Pattern data saved in %s.mat\n',sname);
clear
