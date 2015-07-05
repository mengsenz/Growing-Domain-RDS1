function [ uNew, xNew, tNew ] = domainRemap( u, x, t, tResolution, growthFun )
%DOMAINREMAP remap the location of x from [0,1] back to real location
%(under growth).
%   input:
%           u:              concentration of reactant
%           x:              fixed x positions (nodes)
%           t:              time vector
%           tResolution:    how many total time points wanted
%           growthFun:      a growth function only depends on x and t
%   output:
%           uNew:           down sampled u
%           xNew:           grown x positions
%           tNew:           down sampled t
%{
created: 2015-06-22, MZ
%}

nt = length(t);
tNewIdx = floor(linspace(1,nt,tResolution));
tNew = t(tNewIdx);
uNew = u(:,tNewIdx);
nx = length(x);
xNew = zeros(nx,tResolution);
for i = 1:tResolution
    xNew(:,i) = x*growthFun(tNew(i));
end

end

