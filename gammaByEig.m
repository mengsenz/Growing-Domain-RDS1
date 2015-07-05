function [ gamma ] = gammaByEig( nHalfWaves, l, a, b, d )
%GAMMABYEIG Summary of this function goes here
%   Detailed explanation goes here

k=nHalfWaves*pi/l;
term1 = ((d-1)/(d+1))^2;
term2 = -1+2*b/(a+b)+(a+b)^2;

c1 = term1*term2^2+1-2*b/(a+b)+5*(a+b)^2;
c2 = 2*k^2*(d-1)*term2*(term1-1);
c3 = (term1*4*(d-1)^2-1+3*d)*k^4;

gamma = roots([c1 c2 c3]);
end

