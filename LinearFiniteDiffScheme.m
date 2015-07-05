function [u,v,epsilonu,epsilonv] = LinearFiniteDiffScheme(a,m,dt,ts,numnode,numelem,W,Qx,element,node,d1,d2,aa,bb,gamma,growthFun)
           
 
u = zeros(numnode,ts);
v = zeros(numnode,ts);
% --- random initial conditions
u(1:numnode,1) = aa + bb + .1*rand(numnode,1); 
v(1:numnode,1) = bb/(aa + bb) + .1*rand(numnode,1);

epsilonu=zeros(1,ts);
epsilonv=zeros(1,ts);
% --- intial RHS
fuv1 = assemble_fuv(numnode,numelem,W,Qx,element,node,u(:,1),v(:,1),aa,gamma);
guv1 =assemble_guv(numnode,numelem,W,Qx,element,node,bb,u(:,1),v(:,1),gamma);

u(:,2) = u(:,1) + dt*fuv1 - dt*d1*a*u(:,1);
v(:,2) = v(:,1) - dt*guv1 - dt*d2*a*v(:,1);
 
     

 
for k=2:ts-1
    k
    % --- u(k-1), v(k-1)
    uk = u(:,k-1);
    vk = v(:,k-1);
    % --- u(k), v(k)
    ukk = u(:,k);
    vkk = v(:,k);
    
    % --- intergration for f(u)*phi
    fuvk = assemble_fuv(numnode,numelem,W,Qx,element,node,uk,vk,aa,gamma);
    fuvkk = assemble_fuv(numnode,numelem,W,Qx,element,node,ukk,vkk,aa,gamma);
    % --- intergration for f(v)*phi
    guvk =assemble_guv(numnode,numelem,W,Qx,element,node,bb,uk,vk,gamma);
    guvkk =assemble_guv(numnode,numelem,W,Qx,element,node,bb,ukk,vkk,gamma);
    
    
    scaling = growthFun(dt*k)^2; % domain scaling parameter (Crampin 1999)
    xu = m + 0.5*dt*d1*a/scaling; 
    xv = m + 0.5*dt*d2*a/scaling;

    u(:,k+1) = xu\(1.5*dt*fuvkk - 0.5*dt*fuvk - 0.5*dt*d1*a*ukk/scaling + m*ukk);
    v(:,k+1) = xv\(1.5*dt*guvkk - 0.5*dt*guvk - 0.5*dt*d2*a*vkk/scaling + m*vkk);

%     epsilonu(k) = norm(u(:,k+1)-u(:,k)); epsilonv(k) = norm(v(:,k+1)-v(:,k));
    %epsilon(k) = max(epsilon1,epsilon2);
    epsilon1 = norm(u(:,k+1)-u(:,k))
    %     epsilon2 = norm(v(:,k+1)-v(:,k))
    if isnan(epsilon1), return; end
end
end