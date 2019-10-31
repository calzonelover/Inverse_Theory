% Straight-ray tracing
%
% Copyright 2012,2019 Chaiwoot Boonyasiriwat. All rights reserved.
%
% Input:
%           P1 = [x1,z1] is the starting point
%           P2 = [x2,z2] is the ending point
%           dx           is grid spacing
%           nz           is the model dimension in the z direction
%
% Output:
%           (x,z) = point intersected the horizontal or vertical grid line
%           l     = distance from point (x,z) to the starting point (x1,z1)

function [l,indx] = ray2(P1,P2,dx,nz)
x1 = P1(1);
z1 = P1(2);
x2 = P2(1);
z2 = P2(2);
xmin = min(x1,x2);
xmax = max(x1,x2);
zmin = min(z1,z2);
zmax = max(z1,z2);
diffx = abs(x1-x2);
diffz = abs(z1-z2);
ix1 = floor(x1/dx)+1;
iz1 = floor(z1/dx)+1;
ix2 = floor(x2/dx)+1;
iz2 = floor(z2/dx)+1;

r = sqrt((x1-x2)^2+(z1-z2)^2);
nmax = floor(2*r/dx);
x = zeros(nmax,1);
z = x;
n = 1;
x(1) = x1;
z(1) = z1;

% Case 1: both points are in the same block
if ix1 == ix2 && iz1 == iz2
%     disp('Case 1');
    n = 2;
    x(2) = x2;
    z(2) = z2;
% Case 2: both points are along a horizontal line
elseif diffx > 1e-10 && diffz < 1e-10
%     disp('Case 2');
    for i=min(ix1,ix2)-1:max(ix1,ix2)+1
        xx = (i-1)*dx;
        if xx >= xmin && xx <= xmax
            n = n+1;
            x(n) = xx;
            z(n) = z1;
        end
    end
% Case 3: both points are along a vertical line
elseif diffx < 1e-10 && diffz > 1e-10
%     disp('Case 3');
    for i=min(iz1,iz2)-1:max(iz1,iz2)+1
        zz = (i-1)*dx;
        if zz >= zmin && zz <= zmax
            n = n+1;
            x(n) = x1;
            z(n) = zz;
        end
    end
% Case 4: both points are along a dipping line
else
%     disp('Case 4');
    m = (z2-z1)/(x2-x1);
    a1 = -m;
    c1 = -m*x1+z1;
    for i=min(ix1,ix2)-1:max(ix1,ix2)+1
        xx = (i-1)*dx;
        zz = -(a1*xx-c1);
        if xx >= xmin && xx <= xmax && zz >= zmin && zz <= zmax
            n = n+1;
            x(n) = xx;
            z(n) = zz;
        end
    end
    for i=min(iz1,iz2)-1:max(iz1,iz2)+1
        zz = (i-1)*dx;
        xx = (c1-zz)/a1;
        if xx >= xmin && xx <= xmax && zz >= zmin && zz <= zmax
            n = n+1;
            x(n) = xx;
            z(n) = zz;
        end
    end
end
n = n+1;
x(n) = x2;
z(n) = z2;
x = x(1:n);
z = z(1:n);
l = zeros(n-1,1);
indx = zeros(n-1,1);
for i=1:n-1
    dxx = x(i+1)-x(i);
    dzz = z(i+1)-z(i);
    l(i) = sqrt(dxx*dxx+dzz*dzz);
    ixc = floor(0.5*(x(i)+x(i+1))/dx)+1;
    izc = floor(0.5*(z(i)+z(i+1))/dx)+1;
    indx(i) = izc + (ixc-1)*nz;
end
