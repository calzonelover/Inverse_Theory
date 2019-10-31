% Fast sweeping method (Zhao, 2004) applied to a homogeneous velocity model
% with v = 1 m/s and spatial step (dx,dy) = 1 m. The point source is at the
% center of the model.
%
% Input:
%   v = velocity model
%   h = grid spacing (dx = dy)
%   sx,sy = physical source location (meter)
%
% Output:
%   T = first-arrival traveltime
%
% Copyright 2019 Chaiwoot Boonyasiriwat. All rights reserved.

function T = fsm(v,h,sx,sy)
[ny,nx] = size(v);
isx = floor(sx/h)+1;
isy = floor(sy/h)+1;
T = 1e6*ones(ny,nx);    % initial time
T(isy,isx) = 0;         % initial time at source point
s = 1./v;               % reciprocal of velocity
sh = s*h;
sh2 = sh.*sh;
threshold = 1e-15;       % threshold value for terminating the sweep
nsweep = 100;
Tnorm = zeros(nsweep,1);
for i=1:nsweep
    % Sweep 1
    for iy=1:ny
        for ix=1:nx
            if ix>1 && ix<nx && iy>1 && iy<ny                % handle inner points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy>1 && iy<ny      % handle left-boundary points
                Tx_min = T(iy,ix+1);
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy>1 && iy<ny     % handle right-boundary points
                Tx_min = T(iy,ix-1);
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif iy == 1 && ix>1 && ix<nx      % handle bottom-boundary points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = T(iy+1,ix);
                dT = Tx_min-Ty_min;
            elseif iy == ny && ix>1 && ix<nx     % handle top-boundary points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = T(iy-1,ix);
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy == 1           % handle lower-left corner
                Tx_min = T(1,2);
                Ty_min = T(2,1);
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy == 1          % handle lower-right corner
                Tx_min = T(1,nx-1);
                Ty_min = T(2,nx);
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy == ny          % handle upper-left corner
                Tx_min = T(ny,2);
                Ty_min = T(ny-1,1);
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy == ny         % handle upper-right corner
                Tx_min = T(ny,nx-1);
                Ty_min = T(ny-1,nx);
                dT = Tx_min-Ty_min;
            end
            if abs(dT) > sh(iy,ix)
                Tnew = min(Tx_min,Ty_min)+sh(iy,ix);
            else
                Tnew = 0.5*(Tx_min+Ty_min+sqrt(2*sh2(iy,ix)-dT*dT));
            end
            T(iy,ix) = min(T(iy,ix),Tnew);
        end
    end
    
    % Sweep 2
    for iy=1:ny
        for ix=nx:-1:1
            if ix>1 && ix<nx && iy>1 && iy<ny                % handle inner points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy>1 && iy<ny      % handle left-boundary points
                Tx_min = T(iy,ix+1);
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy>1 && iy<ny     % handle right-boundary points
                Tx_min = T(iy,ix-1);
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif iy == 1 && ix>1 && ix<nx      % handle bottom-boundary points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = T(iy+1,ix);
                dT = Tx_min-Ty_min;
            elseif iy == ny && ix>1 && ix<nx     % handle top-boundary points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = T(iy-1,ix);
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy == 1           % handle lower-left corner
                Tx_min = T(1,2);
                Ty_min = T(2,1);
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy == 1          % handle lower-right corner
                Tx_min = T(1,nx-1);
                Ty_min = T(2,nx);
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy == ny          % handle upper-left corner
                Tx_min = T(ny,2);
                Ty_min = T(ny-1,1);
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy == ny         % handle upper-right corner
                Tx_min = T(ny,nx-1);
                Ty_min = T(ny-1,nx);
                dT = Tx_min-Ty_min;
            end
            if abs(dT) > sh(iy,ix)
                Tnew = min(Tx_min,Ty_min)+sh(iy,ix);
            else
                Tnew = 0.5*(Tx_min+Ty_min+sqrt(2*sh2(iy,ix)-dT*dT));
            end
            T(iy,ix) = min(T(iy,ix),Tnew);
        end
    end

    % Sweep 3
    for iy=ny:-1:1
        for ix=nx:-1:1
            if ix>1 && ix<nx && iy>1 && iy<ny                % handle inner points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy>1 && iy<ny      % handle left-boundary points
                Tx_min = T(iy,ix+1);
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy>1 && iy<ny     % handle right-boundary points
                Tx_min = T(iy,ix-1);
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif iy == 1 && ix>1 && ix<nx      % handle bottom-boundary points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = T(iy+1,ix);
                dT = Tx_min-Ty_min;
            elseif iy == ny && ix>1 && ix<nx     % handle top-boundary points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = T(iy-1,ix);
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy == 1           % handle lower-left corner
                Tx_min = T(1,2);
                Ty_min = T(2,1);
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy == 1          % handle lower-right corner
                Tx_min = T(1,nx-1);
                Ty_min = T(2,nx);
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy == ny          % handle upper-left corner
                Tx_min = T(ny,2);
                Ty_min = T(ny-1,1);
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy == ny         % handle upper-right corner
                Tx_min = T(ny,nx-1);
                Ty_min = T(ny-1,nx);
                dT = Tx_min-Ty_min;
            end
            if abs(dT) > sh(iy,ix)
                Tnew = min(Tx_min,Ty_min)+sh(iy,ix);
            else
                Tnew = 0.5*(Tx_min+Ty_min+sqrt(2*sh2(iy,ix)-dT*dT));
            end
            T(iy,ix) = min(T(iy,ix),Tnew);
        end
    end

    % Sweep 4
    for iy=ny:-1:1
        for ix=1:nx
            if ix>1 && ix<nx && iy>1 && iy<ny                % handle inner points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy>1 && iy<ny      % handle left-boundary points
                Tx_min = T(iy,ix+1);
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy>1 && iy<ny     % handle right-boundary points
                Tx_min = T(iy,ix-1);
                Ty_min = min(T(iy-1,ix),T(iy+1,ix));
                dT = Tx_min-Ty_min;
            elseif iy == 1 && ix>1 && ix<nx      % handle bottom-boundary points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = T(iy+1,ix);
                dT = Tx_min-Ty_min;
            elseif iy == ny && ix>1 && ix<nx     % handle top-boundary points
                Tx_min = min(T(iy,ix-1),T(iy,ix+1));
                Ty_min = T(iy-1,ix);
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy == 1           % handle lower-left corner
                Tx_min = T(1,2);
                Ty_min = T(2,1);
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy == 1          % handle lower-right corner
                Tx_min = T(1,nx-1);
                Ty_min = T(2,nx);
                dT = Tx_min-Ty_min;
            elseif ix == 1 && iy == ny          % handle upper-left corner
                Tx_min = T(ny,2);
                Ty_min = T(ny-1,1);
                dT = Tx_min-Ty_min;
            elseif ix == nx && iy == ny         % handle upper-right corner
                Tx_min = T(ny,nx-1);
                Ty_min = T(ny-1,nx);
                dT = Tx_min-Ty_min;
            end
            if abs(dT) > sh(iy,ix)
                Tnew = min(Tx_min,Ty_min)+sh(iy,ix);
            else
                Tnew = 0.5*(Tx_min+Ty_min+sqrt(2*sh2(iy,ix)-dT*dT));
            end
            T(iy,ix) = min(T(iy,ix),Tnew);
        end
    end
    Tnorm(i) = norm(T(:));
    if i > 1
        if abs(Tnorm(i)-Tnorm(i-1)) < threshold     % termination condition
            break;
        end
    end
end
% disp(['Use ' num2str(i) ' sweeps']);
