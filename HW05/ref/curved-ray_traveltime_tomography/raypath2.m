% Raytracing based on traveltime map
%
% Inputs:
%   t = traveltime table
%   dx = grid spacing
%   s = source position (sx,sz)
%   g = receiver position (gx,gz)
%
% Copyright 2011 Chaiwoot Boonyasiriwat. All rights reserved.
%
% Written by Chaiwoot Boonyasiriwat (July 27, 2011)
% Last modified on Oct 16, 2019

function L = raypath2(t,dx,s,g)
debug = 0;
[nz,nx] = size(t);
xs = s(1);
zs = s(2);

% Raytracing from receiver back to source
x = g(1);
z = g(2);
r = sqrt((x-xs)^2+(z-zs)^2);
L = zeros(1,nx*nz);

if debug
    imagesc((0:nx-1)*dx,(0:nz-1)*dx,t);
    hold on;plot(xs,zs,'wp',x,z,'wv');
end

while r > 1.5*dx
    ix = floor(x/dx)+1;
    iz = floor(z/dx)+1;
    if ix < 2 || ix > nx-1 || iz < 2 || iz > nz-1
        return;
    end
    
    % Compute gradient at x
    [gx,gz] = grad_t(t(iz-1:iz+1,ix-1:ix+1));
    if isnan(gx) || isnan(gz)
        disp('grad is NaN');
        gx = x-xs;
        gz = z-zs;
        g = sqrt(gx*gx+gz*gz);
        gx = gx/g;
        gz = gz/g;
    end
    
    % Move to a new location
    [l,indx] = ray2([x z],[x-gx*dx, z-gz*dx],dx,nz);
    if ~isempty(l)
        for i=1:length(l)
            L(indx(i)) = L(indx(i)) + l(i);
        end
    end
    
    x = x - gx*dx;
    z = z - gz*dx;
    r = sqrt((x-xs)^2+(z-zs)^2);
    
    if debug
%         plot(x,z,'w*');drawnow;
    end
end

% Move to a new location
[l,indx] = ray2([x z],[xs zs],dx,nz);
if ~isempty(l)
    for i=1:length(l)
        L(indx(i)) = L(indx(i)) + l(i);
    end
end

sqrt2 = sqrt(2);
for i=1:length(L)
    if L(i) > dx*sqrt2
        L(i) = dx*sqrt2;
%         error(['too long:' num2str(dx*sqrt2) ',' num2str(L(i))]);
    end
end

if debug
    hold off
    imagesc(reshape(L,nz,nx));drawnow;
end
