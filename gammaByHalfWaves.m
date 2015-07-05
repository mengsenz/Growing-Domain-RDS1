function [ gamma ] = gammaByHalfWaves( nWanted, a,b,d,l, kineticModel, kGM)
%GAMMAHALFPERIOD Summary of this function goes here
%   Detailed explanation goes here
if nargin<7,kGM = 0; end % Gierer-Meinhardt activator saturation
if nargin<6, kineticModel =1; end
switch kineticModel
    case 1 % Schnakenberg
        gamma = 2*nWanted^2*pi^2/l^2/(2*b/(a+b)-1-(a+b)^2/d);
    case 2 % Gierer-Meinhardt
        [ueq, veq] = GiererMeinhardtEquilibrium(a,b,kGM);
        fu = -b+2*ueq/(veq*(1+kGM*ueq^2))-2*kGM*ueq^3/(1+kGM*ueq^2)/veq;
%         fv = -ueq^2/(1+kGM*ueq^2)/veq^2;
%         gu = 2*ueq;
        gv = -1;
        gamma = 2*d/(d*fu+gv)*nWanted^2*pi^2/l^2;
end
        

end

