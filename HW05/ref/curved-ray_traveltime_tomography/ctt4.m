% Driver program for curved-ray traveltime tomography
%
% Copyright 2019 Chaiwoot Boonyasiriwat. All rights reserved.

clear;clc;
itermax = 10;   % maximum number of iteration
nx = 461;
nz = 151;
fid = fopen('vel_nx461_nz151_dx20.dat','r');
v = reshape(fread(fid,nx*nz,'float32'),nz,nx);
dx = 20;
x = (0:nx-1)*dx;
z = (0:nz-1)*dx;
xmax = max(x);
zmax = max(z);
s = 1./v;  % slowness

% sources and receivers
sx = 2*dx:10*dx:xmax-2*dx;
ns = length(sx);
sz = 2*dx*ones(ns,1);
rx = 2*dx:10*dx:xmax-2*dx;
nr = length(rx);
rz = 3*dx*ones(nr,1);
L = raymatrix2(v,dx,sx,sz,rx,rz);
t_true = L*s(:);
n = 1e-2*randn(size(t_true));
delta = norm(n);  % noise level
t_obs = t_true+n;
subplot(231);plot(t_obs,'b');
xlabel('data trace');ylabel('traveltime (s)');
title('Traveltime data');
subplot(235);imagesc(v,[1500 5500]);
title('True Velocity');
colormap(jet(256));
v0 = unpadmodel(diffuse2d(padmodel(v,20),10000,500),20);
subplot(232);imagesc(v0,[1500 5500]);title('Initial Velocity');
s0 = 1./v0;
residual = zeros(itermax,1);
tic
for iter=1:itermax
    L = raymatrix2(v0,dx,sx,sz,rx,rz);
    r = L*s0(:) - t_obs;
    res = norm(r);
    residual(iter) = res;
    disp(['res = ' num2str(res)]);
    g = reshape(L'*r,nz,nx);
    g = unpadmodel(diffuse2d(padmodel(g,20),100,100),20);
    g = g/max(abs(g(:)));
    subplot(233);imagesc(g,[-1,1]);title('Gradient');colorbar
    % back-tracking line search
    step = 1e-4;
    for i=1:20
        s1 = s0 - step*g;
        v1 = 1./s1;
        L = raymatrix2(v1,dx,sx,sz,rx,rz);
        r = L*s1(:) - t_obs;
        res_trial = norm(r);
        disp(['res_trial = ' num2str(res_trial) ', step = ' num2str(step)]);
        if res_trial < res
            break;
        else
            step = 0.5*step;
        end
    end
    s0 = s1;
    v0 = v1;
    subplot(234);plot(residual(1:iter));title('Residual');
    subplot(236);imagesc(v0,[1500 5500]);
    title('Reconstructed Velocity');drawnow
end
residual(itermax+1) = res_trial;
toc