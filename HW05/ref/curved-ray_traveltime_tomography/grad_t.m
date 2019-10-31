% Computing a gradient of 3x3 traveltime matrix using central finite difference
%
% Copyright 2019 Chaiwoot Boonyasiriwat. All rights reserved.

function [gx,gz]=grad_t(t)
gx = t(2,3)-t(2,1);
gz = t(3,2)-t(1,2);
sizeg = sqrt(gx*gx+gz*gz);
gx = gx/sizeg;
gz = gz/sizeg;
