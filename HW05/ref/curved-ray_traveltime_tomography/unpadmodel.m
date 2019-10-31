% 2D unpadding of a model
%
% Copyright 2019 Chaiwoot Boonyasiriwat. All rights reserved.

function v = unpadmodel(v,w)
v = v(w+1:end-w,w+1:end-w);