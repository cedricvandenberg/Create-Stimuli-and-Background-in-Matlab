function GetBackground;
%generate Dirichlet distributions with clumped foci so size is variable
%utilize relative areas of 13 color classes
clear;
%Read 33 background color RGB triplets and proportions
fid=fopen('CedricsRGBproportions.csv');
rd=textscan(fid,'%s %s %s %s %s %s %s %s',1,'delimiter',','); heads=[rd{1:end}];
rd=textscan(fid,'%s %s %f %f %f %f %f %f','delimiter',',');
fclose(fid);
groupname=rd{2}; clrs=[rd{3:5}]; group=rd{6}; lggroup=rd{7}; gprop=rd{8};
clear fid rd; n=length(groupname);

%get proportions of the 13 groups
ngr=35;
relarea=zeros(ngr,1);
for g=1:ngr
  relarea(g)=sum(gprop(group==g));
end;
%get cumulaive probabilities
cumprob=zeros(ngr,1); cumprob(1)=relarea(1); tot=cumprob(1);
for k=2:ngr
  tot=tot+relarea(k);
  cumprob(k)=tot;
end;
clear k tot;

fprintf(1,'Read background RGB triplets, creating mosaic (about 3 min)\n');
tic

np=15000;   %number of patches, 15000 works too

%random foci tesselation
crad=0.3*sqrt(1/np);
%first generate n/4 random points
% n2=round(np/4);
% n2=round(np/5);
n2=round(np/6);
mx=rand(n2,2); x=mx(:,1); y=mx(:,2); clear mx;
%now generate the rest which are closer to existing points than crad
i=n2;
while i<=np
  nx=rand; ny=rand; dst=ones(n2,1);
  for t=1:i
    tx=x(t); ty=y(t);
    dst(t)=sqrt(((nx-tx)^2)+((ny-ty)^2));
  end; %t
  if min(dst)<crad
      i=i+1;
      x(i)=nx; y(i)=ny;
  end;
end;
%figure; plot(x,y,'.b'); axis equal; title([num2str(np) ' points, crad coeff=0.3']);
clear crad n2 mx nx ny i t tx ty dst;


% figure; set(gcf,'Position',[660 135 1256 997]); hold on;
%patch([0 1 1 0],[0 0 1 1],[0.7 0.7 0.7],'EdgeColor','none');
X=[x y]; 
[V,C]=voronoin(X); [nv,i]=size(C); xx=[]; yy=[]; nt=0;
for i=1:nv
  js=C{i}; [j,nj]=size(js); tx=zeros(nj,1); ty=tx;
  for j=1:nj
    tx(j)=V(js(j),1); ty(j)=V(js(j),2);
  end;
  mnx=min(tx); mxx=max(tx); mny=min(ty); mxy=max(ty); 
  if mnx>=0 && mxx<=1 && mny>=0 && mxy<=1
    nt=nt+1; tx=[tx; tx(1)]; ty=[ty; ty(1)]; %close ends
    [a,b,c,d,cnx,cny]=EFourier(tx,ty,10);    %round corners
    [tx,ty,ncu]=InvEFourier(a,b,c,d,cnx,cny,3); clear a b c d cnx cny ncu;
    xx{nt}=tx; yy{nt}=ty;
    % clr=[rand(1,1) rand(1,1) rand(1,1)]; 
    % c=rand(1,1); clr=[c c c]; %to make grays
    % clrs{nt}=clr;
%     plot(tx,ty,'k'); %does outline
%    patch(tx,ty,clr,'EdgeColor','none');  %does a colored patch within the outline
  end;
  clear js j mnx mxx mny mxy nj tx ty; 
end; %vertex 1 to nv
clear V C i x y X i nv; 
%nt is the number of tesselations (smaller than np), each one has 200 points
% axis equal; axis off;
% title([num2str(np) ' clumped foci']); drawnow;

%all patches are stored in xx{nt},yy{nt} where nt = number of tesselations (patches)

%find centroids of each patch
cx=zeros(nt,1); cy=cx;
for t=1:nt
    px=xx{t}; py=yy{t}; [pcx,pcy]=polycenter(px,py);
    cx(t)=pcx; cy(t)=pcy;
end;

% save data
% halt; % **************************************************************************

% figure; set(gcf,'Position',[690 32 1226 1098]); 
% plot(cx,cy,'.k'); hold on; %centroids
% for t=1:nt
%    px=xx{t}; py=yy{t}; plot(px,py,'k');
% end;
% axis equal; axis off;
% title([num2str(nt) ' patches']); clear px py t;

ng=1000;
%choose ng patches at random as nuclei for aggregations of colours in same lggroup 
% and assign 1-13 color classes at random according to relarea-->cumprob

knn=10; %number of nearest neighbors 24 10

aidx=round(nt*rand(ng,1)); aidx(aidx<1)=1; aidx(aidx>nt)=nt;
fx=cx(aidx); fy=cy(aidx);
% plot(fx,fy,'.b','MarkerSize',20); %foci to start group
edges=[0; cumprob]; rs=rand(ng,1); fclr=zeros(ng,1);
for g=1:ng
  hc=histcounts(rs(g),edges); idx=find(hc>0);
  fclr(g)=idx;
end;
clear rs h hc idx; % fclr(ng) contains color codes for each focal patch

patchcolor=zeros(nt,3); %patch color triplet
done=zeros(nt,1); % 1 when done to avoid overlaps
for p=1:ng  
  tx=fx(p); ty=fy(p); pidx=aidx(p); %focus of patch, x,y, index in full list
  [nx,ny]=Circle(tx,ty,0.2);
  in=inpolygon(cx,cy,nx,ny); %1= in the circle, 0 not, size of in =nt
  idx=find(in); %indices to all foci in the circle including patch focus tx,ty
  ngx=cx(idx); ngy=cy(idx); %coordinates of points inside the circle including focus
  %find row in idx with (idx=pidx) with patch focus
  ridx=find(idx==pidx);
  ngbrs=[ngx ngy]; D=pdist(ngbrs); Z=squareform(D);
  row=Z(ridx,:); s=length(row); row=[row' (1:s)'];
  %get the 5 nearest neighbors
  srt=sortrows(row,1); idxs=srt(2:end,2); nn=idxs(1:knn);
  nnidx=idx(nn); nnx=cx(nnidx); nny=cy(nnidx); 
  %nnx and nny are the knn nearest foci
%   figure; set(gcf,'Position',[1066 430 850 700]);
%   plot(tx,ty,'.b','MarkerSize',15);
%   hold on; plot(nx,ny,'b');
%   axis equal
%   plot(cx(cx>(tx-0.055) & cx<(tx+0.055) & cy>(ty-0.055) & cy<(ty+0.055)),...
%        cy(cx>(tx-0.055) & cx<(tx+0.055) & cy>(ty-0.055) & cy<(ty+0.055)),'k.');
%   plot(ngx,ngy,'ok'); plot(nnx,nny,'sb','MarkerSize',8);
%   clear nx ny in idx ngx ngy rids ngbrs D Z row s srt idxs nn;
%   x=xx{pidx}; y=yy{pidx}; plot(x,y,'k');
%   for k=1:knn
%     x=xx{nnidx(k)}; y=yy{nnidx(k)}; plot(x,y,'g');
%   end;
  %VARIABLES FOR THIS FOCUS PATCH ***************************************************
  % xx,yy{1:nt}        patch outlines for all nt patches and centroids
  % tx,ty,pidx         focus patch centroid and index to all centroids and patches  
  % nnx,nny,nnidx(knn) knn nearest patch centroids and index

%   % show patches
%   x=xx{pidx}; y=yy{pidx}; patch(x,y,'b','FaceAlpha',0.5);
%   for k=1:knn
%     xn=xx{nnidx(k)}; yn=yy{nnidx(k)}; patch(xn,yn,'g','Facealpha',0.5);
%   end;
%   
  %get colours for current patch group using clrs(50,3) (rgb per row, 50 rows) 
  %   group(50) (13 color groups) 
  cg=fclr(p); %color group (1:13) for this focal patch.  Get all RGBs for this group
  rgb=clrs(group==cg,:); [npcl,c]=size(rgb); %number of colors in color group cg 
  cclst=round(npcl*rand((knn+1),1)); cclst(cclst==0)=1; cclst(cclst>npcl)=npcl;
  %cclst is index to rows of rgb; size knn+1 (last one is focal color)
  pltclrs=rgb(cclst,:);
  patchcolor(pidx,:)=pltclrs(knn+1,:); 
  done(pidx)=1; done(nnidx)=1;
  for k=1:knn
    patchcolor(nnidx(k),:)=pltclrs(k,:);
  end;
end; %  focal patch and patch group p

%now give random colors to remaining patches (done==0) proportional to areas
ndone=find(sum(patchcolor,2)==0); nnd=length(ndone); % indices to not done patches
rs=rand(nnd,1); %nnd random numbers
edges=[0; cumprob];  rfclr=zeros(nnd,1);
for k=1:nnd
  hc=histcounts(rs(k),edges); idx=find(hc>0);
  rfclr(k)=idx;
end;
clear k rs hc idx; % refclr(nnd) has color codes 1-13 for each of nnd single patches
%collect RGBs for each color group
for k=1:nnd
  cg=rfclr(k); 
  rgb=clrs(group==cg,:); [npcl,c]=size(rgb); %number of colors in color group cg 
  if npcl>1
      cc=round(5*rand); cc(cc==0)=1; cc(cc>npcl)=npcl;
      %cc is randomly drawn color code in rgb for this group
  else
      cc=1;
  end;
  patchcolor(ndone(k),:)=rgb(cc,:);
end;
ttl=[num2str(nt) ' patches, ' num2str(ng) ' groups, ' num2str(knn) ' neighbors/group'];
clear D Z aidx c* done e* f* g* i* k* lggroup ndone ng ngbrs ngx ngy nn nnd nnidx;
clear nx ny nnx nny np npcl ox oy p pcx pcy pidx pltclrs px py r* s* t tx ty x y;

%save the background as a .mat file
prompt={'Enter background data file name to save'}; tt='BACKGROUND FILE NAME';
fn=inputdlg(prompt,tt,[1 30]); fn=char(fn); fn=['BACK' fn];
save(fn);
fprintf(1,'Background saved in %s\n',fn);

fprintf(1,'Making background image display, takes a minute to show\n');
figure; set(gcf,'Position',[90 32 1226 1098]);
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
fprintf(1,'Finished background\n');
end

function [a,b,c,d,cnx,cny]=EFourier(x,y,nc)
% [a,b,c,d,cnx,cny]=EFourier(x,y,nc);
% Function to calculate elliptical fourier coefficients for an
%   outline.  Outline as vectors of (x,y) points.  
% x,y  Input vectors, each row is a point
% nc   Number of coefficients to calculate (better fit with more)
%         If none specified then will take nc=40
% a,b,c,d  Elliptical Fourier Coefficients (real and imaginary parts)
% cnx,cny  Geometrical center in original coordinates
% use InvEFourier to recover shape from the coefficients
% for PCA use  c1=complex(a,b); c2=complex(c,d); c1=abs(c1); c2=abs(c2);           
% John A. Endler, 2004, corrections: John.Endler@jcu.edu.au
%
tpi=2*pi;
if nargin==2 nc=40; end;
[n,i]=size(x); if n==1 x=x'; y=y'; end; %convert to column vectors first
%remove duplicate points but keep point order
dta=[x y]; dta=unique(dta,'rows','stable'); x=dta(:,1); y=dta(:,2); clear dta;
[n,i]=size(x); %number of points
%find geometrical center and move object so that this center is at (0,0)
[cx,cy,cnx,cny]=GeoCenter(x,y);
x=cx; y=cy;

%convert (x,y) to differences between points
for i=2:n dx(i-1)=x(i)-x(i-1); dy(i-1)=y(i)-y(i-1); end;
dx(n)=x(1)-x(n); dy(n)=y(1)-y(n);
tsum=0; xi=0; yi=0; delta=0; xsum=0; ysum=0; ao=0; co=0;
for i=1:n
  t(i)=sqrt((dx(i)^2)+(dy(i)^2));
  if t(i)>0
    rdx(i)=dx(i)/t(i);
    rdy(i)=dy(i)/t(i);
    tnew=tsum+t(i);
    xi=xsum-rdx(i)*tsum;
    yi=ysum-rdy(i)*tsum;
    t1=tnew-tsum;
    t2=(tnew^2)+(tsum^2);
    ao=ao+(rdx(i)*t2/2)+xi*t1;
    co=co+(rdy(i)*t2/2)+yi*t1;
    tsum=tnew;
  else
    rdx(i)=0; rdy(i)=0; 
  end;
end;
ao=ao/tsum; co=co/tsum; tlen=tsum;

for h=1:nc
  fh=h;
  fact1=fh*tpi/tsum;
  angprv=0; asum=0; bsum=0; csum=0; dsum=0;
  for i=1:n
     ang=angprv+t(i)*fact1;
     wtac=cos(ang)-cos(angprv);
     wtbd=sin(ang)-sin(angprv);
     angprv=ang;
     asum=asum+rdx(i)*wtac;
     bsum=bsum+rdx(i)*wtbd;
     csum=csum+rdy(i)*wtac;
     dsum=dsum+rdy(i)*wtbd;
  end;
  fact2=fact1*fh*pi;
  a(h)=asum/fact2;
  b(h)=bsum/fact2;
  c(h)=csum/fact2;
  d(h)=dsum/fact2;
end;
a=a'; b=b'; c=c'; d=d';  %transpose to column vectors
end %efourier

function [px,py,ncu]=InvEFourier(a,b,c,d,cnx,cny,nc,nang)
% [px,py,ncu]=InvEFourier(a,b,c,d,cnx,cny,nc,nang);
% Recovers shape outline from fourier coefficients
% INPUT
%   a b c d   fourier coefficient column vectors from EFourier
%   cnx cny   geometrical center of original outline from EFourier
%               Set these to 0 for reconstruction centered at (0,0)
%   nc        number of coefficients to use in reconstruction
%   nang      number of angles to use in reconstructin (default 200)
% OUTPUT
%   px py     column vectors with the reconstructed outline
%   ncu       number of coefficients used (<=size of a b c & d)
% John A. Endler, 2004, corrections: John.Endler@jcu.edu.au
%
if nargin<8 n=200; else n=nang; end;
mnc=size(a); if nc>mnc nc=mnc; end; ncu=nc;
tpi=2*pi;
%number of angles to use around outline
    for i=1:n
      xi=0; yi=0;
      delang=tpi*(i-1)/n; ang=delang;
      for h=1:ncu
        ca=cos(ang); sa=sin(ang);
        xi=xi + a(h)*ca + b(h)*sa;
        yi=yi + c(h)*ca + d(h)*sa;
        ang=ang+delang;
      end;
      rx(i)=xi; ry(i)=yi;
    end;
    rx=rx'; ry=ry';
    [cx,cy,tcx,tcy]=GeoCenter(rx,ry);
    px=cx; py=cy; clear cx cy tcx tcy;
    px=px+cnx; py=py+cny;
end %inv Efourier

function [cx,cy,cnx,cny]=GeoCenter(x,y);
% [cx,cy,cnx,cny]=GeoCenter(x,y);
% Function to calculate geometric center of an outline
%  and move it so that center is at origin (0,0)
% x,y    Input vectors, each row is a point
% cx,cy   Output coordinates, centered at (0,0)
% cnx,cny Geometric center in original coordinates
[n,i]=size(x); %number of points
%first find center
for j=1:n
  i=j+1; if j==n i=1; end;
  aa(j)=((y(i)+y(j))*(x(j)-x(i)))/2;
  my(j)=((y(i)^2)+y(i)*y(j)+(y(j)^2))*(x(j)-x(i))/6;
  mx(j)=((x(i)^2)+x(i)*x(j)+(x(j)^2))*(y(j)-y(i))/6;
end;
sx=0; sy=0; sa=0;
for j=1:n sa=sa+aa(j); sx=sx+mx(j); sy=sy+my(j); end;
cnx=-sx/sa; cny=sy/sa;  %geometric center in original coordinates
% now move such that geometric center is at 0,0
for i=1:n cx(i)=x(i)-cnx; cy(i)=y(i)-cny; end;
cx=cx'; cy=cy';
end %geocenter

function [cx,cy]=Circle(ox,oy,rd)
% [cx,cy]=Circle(ox,oy,rd);
% draw a circle of radius rd centered on (ox,oy)
r=1;
x=-1:0.01:1; x2=fliplr(x);
y1=(r^2-x.^2).^(1/2);
y2=-(r^2-x.^2).^(1/2);
cy=rd*[y1 y2]; cx=rd*[x x2];
cy=cy+oy; cx=cx+ox;
end %circle
