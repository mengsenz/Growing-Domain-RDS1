function [ lt ] = logisticGrowth( tnow, l0, kappa, r)
%LOGISTICGROWTH Summary of this function goes here
%   [ lt ] = LogisticGrowth( tnow, l0, kappa, r)
lt = l0*kappa*exp(r*tnow)/(kappa-1+exp(r*tnow));
end

