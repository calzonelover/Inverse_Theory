% 1D smoothing filter based on the diffusion equation.
%
% Copyright 2019 Chaiwoot Boonyasiriwat. All rights reserved.

function [x,d,u] = diffuse1d(x,alpha)

if nargin == 1
    alpha = 100;
end;

n = length(x);
a = -alpha;
b = 1-2*a;
c = alpha*alpha;

% Form the sparse system of equation
d = b*ones(n,1); % Main diagonal
d(1) = -c;
d(end)= -c;
u = x;

% Forward elimination
fac = a/d(1);
d(2) = d(2) - c*fac;
u(2) = u(2) - u(1)*fac;
for i=3:n-1
    fac = a/d(i-1);
    d(i) = d(i) - a*fac;
    u(i) = u(i) - u(i-1)*fac;
end;
fac = c/d(end-1);
d(end) = d(end) - a*fac;
u(end) = u(end) - u(end-1)*fac;

% Backward substitution
xx = zeros(n,1);
xx(end) = u(end)/d(end);
for i=n-1:-1:2
    xx(i) = (u(i) - a*xx(i+1))/d(i);
end;
xx(1) = (u(1) - c*xx(2))/d(1);
x = xx;
