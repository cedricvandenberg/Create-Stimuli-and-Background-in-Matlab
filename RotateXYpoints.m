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
