function [ u_eq, v_eq ] = GiererMeinhardtEquilibrium( a, b, k )
%GIERERMEINHARDTEQUILIBRIUM Equilibrium solution without diffusion term for
%Gierer-Meinhardt Model
%   [ u_eq, v_eq ] = GiererMeinhardtEquilibrium( a, b, k )

u_eq = roots([k*b -k*a b -a-1]);
u_eq = u_eq(u_eq>0); % only keep the real positive root
v_eq = u_eq^2;
end

