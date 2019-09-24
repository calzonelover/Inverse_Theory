clear;clc
show_ray = 0;
nx = 50;
nz = 50;
fid = fopen('vel_nx50_nz50_dx20.dat','r');
v = reshape(fread(fid,nx*nz,'float32'),nz,nx);
dx = 20;
x = (0:nx-1)*dx;
z = (0:nz-1)*dx;
xmax = max(x);
zmax = max(z);
s = 1./v;  % slowness

% sources and receivers
sz = 0:dx:zmax;
ns = length(sz);
sx = zeros(ns,1);
rz = sz;
nr = length(rz);
rx = xmax*ones(nr,1);
L = zeros(ns*nr,nx*nz);

for is=1:ns
    for ir=1:nr
        [l,indx,xp,zp,r] = ray(nx+1,nz+1,dx,sx(is),sz(is),rx(ir),rz(ir));
        if show_ray
            subplot(211);imagesc(x,z,v,[1400 4600]);colorbar;axis image
            colormap(jet)
            hold on;plot(sx(is),sz(is),'rp',rx(ir),rz(ir),'rv','MarkerSize',10);
            plot(xp,zp,'w')
            hold off
        end

        LL = zeros(1,nz*nx);
        for i=1:length(l)
            LL(indx(i)) = l(i);
        end
        if show_ray
            subplot(212);imagesc(reshape(LL,nz,nx));axis image
            drawnow
        end
        L(ir+(is-1)*nr,:) = LL;
    end
end
t = L*s(:);
n = 1e-4*randn(size(t));
delta = norm(n);  % noise level
t = t+n;
subplot(211);plot(t,'b');
set(gca,'FontSize',15);
xlabel('data trace')
ylabel('traveltime (s)');
title('Traveltime data');

LL = L'*L;
alpha = 1;
s2 = (LL+alpha*eye(size(LL)))\(L'*t);
subplot(223);imagesc(x,z,v,[1500 2000]);colorbar;axis image
set(gca,'FontSize',15,'XTick',[0 500 900]);
title('True velocity');
subplot(224);imagesc(x,z,reshape(1./s2,nz,nx),[1500 2000]);axis image;colorbar
set(gca,'FontSize',15,'XTick',[0 500 900]);
title('Inverted velocity: \alpha = 1');
colormap(jet)
