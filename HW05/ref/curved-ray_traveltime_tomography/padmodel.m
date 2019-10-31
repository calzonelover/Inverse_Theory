% 2D padding of a model
%
% Copyright 2019 Chaiwoot Boonyasiriwat. All rights reserved.

function vpad = padmodel(v,npad)

[nz,nx] = size(v);
npad1 = npad+1;
nx2 = nx + 2*npad;
nz2 = nz + 2*npad;
vpad = zeros(nz2,nx2);
vpad(npad1:npad+nz,npad1:npad+nx) = v;
for i=1:npad
    vpad(:,i) = vpad(:,npad1);
    vpad(:,nx2-i+1) = vpad(:,nx2-npad);
end;
for i=1:npad
    vpad(i,:) = vpad(npad1,:);
    vpad(nz2-i+1,:) = vpad(nz2-npad,:);
end;

