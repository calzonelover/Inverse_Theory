% This function computes a sparse matrix containing the raypath length
% in each cell from all pairs of sources and receivers.
%
% Copyright 2019 Chaiwoot Boonyasiriwat. All rights reserved.

function L = raymatrix2(v,dx,sx,sz,rx,rz)
[nz,nx] = size(v);
ns = length(sx);
nr = length(rx);
row = zeros(ns*nr,1);
col = row;
val = row;
n = 0;
for is=1:ns
    T = fsm(v, dx, sx(is), sz(is));
    for ir=1:nr
        LL = raypath2(T,dx,[sx(is),sz(is)],[rx(ir),rz(ir)]);
        i = find(abs(LL)>1e-10);
        if ~isempty(i)
            for j=1:length(i)
                n = n +1;
                row(n) = ir+(is-1)*nr;
                col(n) = i(j);
                val(n) = LL(i(j));
            end
        end
    end
end
if n < ns*nr
    row = row(1:n);
    col = col(1:n);
    val = val(1:n);
end
L = sparse(row,col,val,ns*nr,nx*nz);