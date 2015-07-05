function fg = GiererMeinhardtKinetics( u,v,a,b,k )
%GiererMeinhardtKinetics Gierer-Meinhardt Kinetics with the form of 
%                           f = a-b*u+u^2/v/(1+k*u^2);      (the activator)
%                           g = u^2 - v;                    (the inhibitor)
%   usage: [f,g] = GiererMeinhardtKinetics( u,v,a,b,k )
%   
fg = [a-b*u+u^2/v/(1+k*u^2), u^2 - v];
end

