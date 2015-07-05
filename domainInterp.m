function [ actualX, value ] = domainInterp( concentration, scaledX, t, growthFunction)
%DOMAININTERP remap domain for [0,1] to actual length and interpolate data
%   Detailed explanation goes here
nx = size(concentration,1);
dx=scaledX(2)-scaledX(1);
T=t(end);
finalLength = growthFunction(T);
nrow = finalLength*nx;
ncol = size(t,2);

actualX = 0:dx:finalLength;
value = zeros(nrow, ncol);

for i = 1:ncol
    len = growthFunction(t(i));
    newX = 0:dx:len;
    nNewX = length(newX);
    value(1:nNewX,i) = interp1(linspace(0,len,nx), concentration(:,i), newX); 
end

end

