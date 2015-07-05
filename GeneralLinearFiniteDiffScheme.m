function [u,v,epsilonu,epsilonv] = GeneralLinearFiniteDiffScheme(a,m,dt,ts,numnode,numelem,W,Qx,element,node,d1,d2,u0,v0,reactionFun,growthFun, stochasticity)
if nargin < 17
    stochasticity = 0;
else
    stochasticity = sqrt(stochasticity);
end
 
u = zeros(numnode,ts);
v = zeros(numnode,ts);
% --- if constant ICs --> expand to a column vector
if numel(u0)==1, u0 = u0*ones(numnode,1); end
if numel(v0)==1, v0 = v0*ones(numnode,1); end
% --- random initial conditions
u(1:numnode,1) = u0 + .1*rand(numnode,1); 
v(1:numnode,1) = v0 + .1*rand(numnode,1);

epsilonu=zeros(1,ts);
epsilonv=zeros(1,ts);
% --- intial RHS
[fuv1, guv1] = assemble_reaction(numnode,numelem,W,Qx,element,node,u(:,1),v(:,1),reactionFun);

u(:,2) = u(:,1) + dt*fuv1 - dt*d1*a*u(:,1);
v(:,2) = v(:,1) - dt*guv1 - dt*d2*a*v(:,1);
 
     

[fuvkk, guvkk] = assemble_reaction(numnode,numelem,W,Qx,element,node,u(:,1),v(:,1),reactionFun);
for k=2:ts-1
    k
    % --- u(k-1), v(k-1)
%     uk = u(:,k-1);
%     vk = v(:,k-1);
    % --- u(k), v(k)
    ukk = u(:,k);
    vkk = v(:,k);
    
    % --- intergration for f*phi and g*phi
    fuvk = fuvkk;
    guvk = guvkk;
%     [fuvk, guvk] = assemble_reaction(numnode,numelem,W,Qx,element,node,uk,vk,reactionFun);
    [fuvkk, guvkk] = assemble_reaction(numnode,numelem,W,Qx,element,node,ukk,vkk,reactionFun);
    
    
    scaling = growthFun(dt*k)^2; % domain scaling parameter (Crampin 1999)
    xu = m + 0.5*dt*d1*a/scaling; 
    xv = m + 0.5*dt*d2*a/scaling;

    u(:,k+1) = stochasticity*(rand(numnode,1)-0.5)+xu\(1.5*dt*fuvkk - 0.5*dt*fuvk - 0.5*dt*d1*a*ukk/scaling + m*ukk);
    v(:,k+1) = stochasticity*(rand(numnode,1)-0.5)+xv\(1.5*dt*guvkk - 0.5*dt*guvk - 0.5*dt*d2*a*vkk/scaling + m*vkk);

%     epsilonu(k) = norm(u(:,k+1)-u(:,k)); epsilonv(k) = norm(v(:,k+1)-v(:,k));
    %epsilon(k) = max(epsilon1,epsilon2);
    epsilon1 = norm(u(:,k+1)-u(:,k))
    %     epsilon2 = norm(v(:,k+1)-v(:,k))
    if isnan(epsilon1), return; end
end
end