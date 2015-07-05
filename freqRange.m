function [ f ] = freqRange(d, gamma, a, b, kineticModel, kGM)
%FREQRANGE Frequency Range or dispersion relation for a specific kinetic
%model (cycles per unit distance)
%   L: domain Length
if nargin<6, kGM=0; end % Gierer-Meinhardt saturation
if nargin<5, kineticModel = 1; end 

switch kineticModel
    case 1
        fu=2*b/(a+b)-1;
%         fv=(a+b)^2;
%         gu=-2*b/(a+b);
        gv=-(a+b)^2;
        detA=(a+b)^2;
    case 2
        [ueq, veq] = GiererMeinhardtEquilibrium(a,b,kGM);
        fu = -b+2*ueq/(veq*(1+kGM*ueq^2))-2*kGM*ueq^3/(1+kGM*ueq^2)^2/veq;
        fv = -ueq^2/(1+kGM*ueq^2)/veq^2;
        gu = 2*ueq;
        gv = -1;
        detA = fu*gv-fv*gu;
end
        
if detA > (d*fu+gv)^2/4/d, warnning('not in Turing Space!');end
charPoly=[d, -gamma*(d*fu+gv), gamma^2*detA]; % det(gamma*A+k^2*D)
r = roots(charPoly);
f=sqrt(r)/(2*pi);

end

