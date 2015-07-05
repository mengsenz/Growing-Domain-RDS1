function [ x ] = linearGrowth( t, x0, r )
%LINEARGROWTH linear growth function
%   t:  current time
%   x0: initial domain size
%   r:  growth rate
x = x0+r*t;

end

