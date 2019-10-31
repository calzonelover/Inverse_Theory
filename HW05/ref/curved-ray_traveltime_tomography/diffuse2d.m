% 2D smoothing filter based on the diffusion equation.
%
% Copyright 2019 Chaiwoot Boonyasiriwat. All rights reserved.

function x = diffuse2d(x,rx,rz)

if nargin == 2
    rz = rx;
end;

[nz,nx] = size(x);

for ix=1:nx
    x(:,ix) = diffuse1d(x(:,ix),rz);
end;

for iz=1:nz
    x(iz,:) = (diffuse1d(x(iz,:),rx))';
end;