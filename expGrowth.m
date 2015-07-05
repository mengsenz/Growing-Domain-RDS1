function lt = expGrowth( tnow, l, r )
%EXPGROWTH Summary of this function goes here
%   usage: x = expGrowth( tnow, l, r )
%   input:
%           tnow:   current time
%           l:      initial domain length
%           r:      growth rate
%   output:
%           lt:     current domain length
lt=l*exp(r*tnow);

end

