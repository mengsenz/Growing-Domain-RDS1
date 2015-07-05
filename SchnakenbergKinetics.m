function fg = SchnakenbergKinetics( u,v, a,b )
%SCHNAKENBERGKINETICS Schnakenberg Kinetics for Reaction-Diffusion Systems
%taking the forms:
%                   f=a-u+u^2*v;
%                   g=b-u^2*v;
%   [ f,g ] = SchnakenbergKinetics( a,b )

fg=[a-u+u.^2.*v, b-u.^2.*v];
end

