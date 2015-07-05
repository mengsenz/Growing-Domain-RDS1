function [ u_eq,v_eq ] = SchnakenbergEquilibrium( a,b )
%SCHNAKENBERGEQUILIBRIUM Summary of this function goes here
%   Detailed explanation goes here
u_eq = a+b;
v_eq = b/(a+b)^2;

end

