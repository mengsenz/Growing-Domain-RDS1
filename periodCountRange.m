function [ nRange ] = periodCountRange( L, d, gamma, a, b, kineticModel, kGM  )
%PERIODCOUNTRANGE Summary of this function goes here
%   Detailed explanation goes here
if nargin<7, kGM=0; end % Gierer-Meinhardt saturation
if nargin<6, kineticModel = 1; end 

f = freqRange(d,gamma,a,b, kineticModel, kGM); 
nRange = f*L;
end

