% fem_example.m

% A simple 1D finite element program in Matlab
%
close all
nx = 100;                %input('input number of grid points in the x-direction\n')


nn = nx;                  % number of nodes
ne = nn-1;                % number of elements
h  = 2/ne;                % length of each element
                          % number of unknowns

% ---- define mesh ---- %
l = 6;
node = linspace(0,l,nn)';         % node coordinates
element = [ 1:(nn-1); 2:nn ]';     % element connectivity


% ---- TURING PARAMETERS ---- %
initialNumberOfHalfWaves = 12;
d1 = 1;       %d1 = input('positive diffusion coefficient sor u'); 
d2 = 18; %18;%10;        %d2 = input('positive diffusion coefficient for v');

% ---- kinetics ---- %
% .......................................
% kineticFunName = 'Schnakenberg';
% kineticIndex = 1;
% aa = 0.1; %0.14;%0.1;
% bb =0.9;%0.16;%0.9;
% reactionFun = @(unow,vnow) gamma*SchnakenbergKinetics(unow,vnow,aa,bb);
% [u0,v0] = SchnakenbergEquilibrium(aa,bb);
% gamma = gammaByHalfWaves(initialNumberOfHalfWaves, aa,bb,d2,l,1 ); %28.2;%113;%190; %60;%600;

% .......................................
kineticFunName = 'Gierer-Meinhardt';
kineticIndex = 2;
aa = 0.1;
bb = 1;
k = 0; % no activator saturation
gamma = gammaByHalfWaves(initialNumberOfHalfWaves, aa,bb,d2,l,kineticIndex,k); %28.2;%113;%190; %60;%600;
reactionFun = @(unow,vnow) gamma*GiererMeinhardtKinetics(unow,vnow,aa,bb,k);
[u0,v0] = GiererMeinhardtEquilibrium(aa,bb,k);

% ---- time step configuration ---- %
T = 3; %100; %10;
t_in = 0;       %t_in = input('input initial time\n');
t_final = 1*T;  %t_final = input('input final time\n');
nt = 2000*T;      %nt = input('input number of time steps\n');

dt = (t_final-t_in)/nt;         %timestep
t = t_in:dt:t_final;  
ts = size(t,2);                 %number of time iterations or nt+1

% ---- domain growth configuration ---- %
r = 0.05; % growth rate
% .......... Growth Functions:
% .......................................
growthFunName = 'No Growth';
growthFun = @(tnow) 1; 
% .......................................
% growthFunName = 'Linear Growth';
% growthFun = @(tnow) linearGrowth(tnow,l,r); 
% .......................................
% growthFunName = 'Exponential Growth';
% expGrowthRate(r,l,t_final); % adjust growth rate to match the final length of linear growth % r = 0.03; %
% growthFun = @(tnow) expGrowth(tnow,l,r); 
% .......................................
% growthFunName = 'Logistic Growth';
% kappa = linearGrowth(t_final,l,r)+0.5; %1.5; % logistic growth maximum
% r = logisticGrowthRate(kappa, r, l, t_final); % adjust growth rate to match the final length of linear growth
% growthFun = @(tnow) logisticGrowth(tnow,l,kappa,r); 

% return
% ----- Batch processing ----- %

batchSize = 30;
u_final = zeros(nx,batchSize);
v_final = zeros(nx,batchSize);
% return
for i = 1:batchSize

    % DEFINE ASSEMLE MATRIX AND RHS
    %
    % Here we define the  data structures
    %   
    %   f - is the RHS vector.  
    %   a - is the sparse assemble matrix.



    % QUADRATURA INTEGRATION 
    % First we need to calculate the quadrature points. 
    % We are using 3 point quadrature rule

    [W,Qx] = quadrature3(element,node,ne);

    % SETTING UP THE ASSEMBLE MATRIX AND RHS VECTOR

    [a,m] =assemble(nn,ne,W,Qx,element,node);

    % solve system
    tic
%     [u,v,~,~] = LinearFiniteDiffScheme(a,m,dt,ts,nn,ne,W,Qx,element,node,d1,d2,aa,bb,gamma,growthFun);
    [u,v,~,~] = GeneralLinearFiniteDiffScheme(a,m,dt,ts,nn,ne,W,Qx,element,node,d1,d2,u0,v0,reactionFun,growthFun);

    toc
    
    % --- save final u and v
    u_final(:,i) = u(:,ts);
    v_final(:,i) = v(:,ts);
end
      
%% ----- save data ----- %%
filename =['BatchResult_' datestr(now, 'yyyy-mm-dd_HH-MM-SS')];
save(['data/', growthFunName, '_', filename],'u_final','v_final', 'l','node','t','d1','d2','aa','bb','gamma');%,'kappa'
%% 
% figure
% % --- initial state
% subplot(1,2,1); hold on
% plot(node,u(:,1),'r-o')
% plot(node,v(:,1),'b-*')
% title('initial concentration')
% legend({'u','v'})
% % final state
% subplot(1,2,2); hold on
% plot(node,u(:,ts),'r-o')
% plot(node,v(:,ts),'b-*')
% title('final concentration')
% legend({'u','v'})
% 
% tResolution = 1000;
% [uNew, xNew, tNew] = domainRemap(u,node,t,tResolution, growthFun);
% [vNew, ~, ~] = domainRemap(v,node,t,tResolution, growthFun);
% 
% 
% figure('position',[1 1 776 690],'PaperPositionMode','auto')
% 
% subplot(2,1,1)
% % pcolor(t,node,u); shading interp; colorbar
% pcolor(tNew,xNew,uNew); shading interp; colorbar; colormap(jet)
% title('u','fontsize',14)
% xlabel('time')
% ylabel('x')
% 
% subplot(2,1,2)
% % pcolor(t,node,v); shading interp; colorbar
% pcolor(tNew,xNew,vNew); shading interp; colorbar; colormap(jet)
% title('v','fontsize',14)
% xlabel('time')
% ylabel('x')
% 
% % ---- graphic annotation
% if strcmp(growthFunName,'Logistic Growth')
%     annotKappa = [', kappa = ' num2str(kappa)];
% else
%     annotKappa = '';
% end
% annotation(gcf,'textbox',...
%     [0.631504922644164 0.489533011272141 0.367088607594937 0.0740740740740737],...
%     'String',{[kineticFunName ' + ' growthFunName],...
%               ['d_1 = ', num2str(d1), ' 	d_2 = ',num2str(d2),' 	r = ', num2str(r), annotKappa],...
%               ['a = ', num2str(aa), ' 	b = ', num2str(bb), ' 	\gamma = ', num2str(gamma)] },...
%     'LineStyle','none',...
%     'FontWeight','bold',...
%     'FontSize',10,...
%     'FontName','Helvetica-Narrow',...
%     'FitBoxToText','off');
% 
% % ---- save graphic
% filename =['result_' datestr(now, 'yyyy-mm-dd_HH-MM-SS')];
% print(['figures/' filename],'-dpng','-r0')
% %% --- save data
% save(['data/', growthFunName '_' filename],'u','v','tNew','xNew','uNew','vNew')
%% ===== plot the distribution of spatial frequency ===== %%
l_final = max(growthFun(t_final),l);
freq=spatialFreq(u_final,l_final);

% --- Histogram
frange = freqRange(d2,gamma,aa,bb,kineticIndex,k)
nbins = 30;
ctrs = linspace(min(frange),max(frange),nbins+1); ctrs = ctrs(1:end-1)+(ctrs(2)-ctrs(1))/2;
counts = hist(freq,ctrs);
figure
plot(ctrs,counts/batchSize)
xlim([frange(2) frange(1)])
title('distribution of final spatial frequency')
xlabel('spatial frequency (cycle/unit distance)')
ylabel('percentage of occurrence')
annotKappa=''
annotation(gcf,'textbox',...
    [0.5 0.5 0.37 0.4],...
    'String',{[kineticFunName ' + ' growthFunName],...
              ['l = ', num2str(l)],...
              ['d_1 = ', num2str(d1), ' 	d_2 = ',num2str(d2),' 	r = ', num2str(r), annotKappa],...
              ['a = ', num2str(aa), ' 	b = ', num2str(bb), ' 	\gamma = ', num2str(gamma)] },...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',10,...
    'FontName','Helvetica-Narrow',...
    'FitBoxToText','off');

print(['figures/batchResults/' filename],'-dpng','-r0')
savefig(gcf,['figures/batchResults/' filename])
%% ===== remap domain of u ===== %%
% [newNode, newU] = domainInterp(u,node,t,@(tnow) linearGrowth(tnow,l,r));
% 
% figure
% pcolor(t,newNode,newU); shading interp
% xlabel('time')
% ylabel('x')
% title('concentration of u under domain growth')

%% ===== remap domain of v ===== %%
% [newNode, newU] = domainInterp(v,node,t,@(tnow) linearGrowth(tnow,l,r));
% 
% figure
% pcolor(t,newNode,newU); shading interp
% xlabel('time')
% ylabel('x')
% title('concentration of v under domain growth')