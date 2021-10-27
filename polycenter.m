function [pcx,pcy]=polycenter(px,py);
% [pcx,pcy]=polycenter(px,py); finds the geometric center (pcx,pcy)
%   of a polygon whose points are in (px,py) column vectors
[np,i]=size(px);
for j=1:np
    i=j+1;  if j==np i=1; end;
    a(j)=((py(i)+py(j))*(px(j)-px(i)))/2;
    my(j)=(((py(i)^2)+py(i)*py(j)+(py(j)^2))*(px(j)-px(i)))/6;
    mx(j)=(((px(i)^2)+px(i)*px(j)+(px(j)^2))*(py(j)-py(i)))/6;
end;
sx=0; sy=0; sa=0;
for j=1:np
    sa=sa+a(j);
    sx=sx+mx(j);
    sy=sy+my(j);
end;
pcx=-sx/sa;
pcy=sy/sa;
    