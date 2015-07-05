function [ lt ] = quadraticGrowth( tnow, tend, r, l0 )
%QUADRATICGROWTH Summary of this function goes here
%   Detailed explanation goes here

lt = l0-r*(tnow-tend)*tnow;
end

