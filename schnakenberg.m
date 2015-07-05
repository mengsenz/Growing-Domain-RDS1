% fem_example.m

% A simple 1D finite element program in Matlab
%
nx = 100;                %input('input number of grid points in the x-direction\n')


nn = nx;                  % number of nodes
ne = nn-1;                % number of elements
h  = 2/ne;                % length of each element
                          % number of unknowns

% ---- define mesh ---- %
l = 4.4781;%testX; %1;
node = linspace(0,l,nn)';         % node coordinates
element = [ 1:(nn-1); 2:nn ]';     % element connectivity


%HOPF BIFURCATION PARAMETERS
% d1 = 1;       %d1 = input('positive diffusion coefficient sor u'); 
% 
% d2 = 0.43;        %d2 = input('positive diffusion coefficient for v');
% 
% alpha = 0.3^3;
% beta = 0.3;
% gamma =1;
% aa = (beta - alpha)/2;
% bb = (beta + alpha)/2;

% Nd = 15;
% ds = linspace(10,100, Nd);
% for j = 1:Nd
% ---- TURING PARAMETERS ---- %
initialNumberOfHalfWaves = 2;
d1 = 1;       %d1 = input('positive diffusion coefficient sor u'); 
d2 = 18; %18;%10; d2 = ds(j);  %d2 = input('positive diffusion coefficient for v');
% return

% ---- kinetics ---- %
% .......................................
kineticFunName = 'Schnakenberg';
kineticIndex = 1;
aa = 0.1; %0.14;%0.1;
bb =0.9;%0.16;%0.9;
gamma = gammaByHalfWaves(initialNumberOfHalfWaves, aa,bb,d2,1,1 ); %28.2;%113;%190; %60;%600;
reactionFun = @(unow,vnow) gamma*SchnakenbergKinetics(unow,vnow,aa,bb);
[u0,v0] = SchnakenbergEquilibrium(aa,bb);
% u0=testU;
% v0=testV;
noiseLevel = 0;%0.01;

% .......................................
% kineticFunName = 'Gierer-Meinhardt';
% kineticIndex = 2;
% aa = 0.1;
% bb = 1;
% k = 0; %0.05; % no activator saturation
% gamma = gammaByHalfWaves(initialNumberOfHalfWaves, aa,bb,d2,l,kineticIndex,k); %28.2;%113;%190; %60;%600;
% reactionFun = @(unow,vnow) gamma*GiererMeinhardtKinetics(unow,vnow,aa,bb,k);
% [u0,v0] = GiererMeinhardtEquilibrium(aa,bb,k);
% noiseLevel = 0; %0.01; %0.05; 


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
% % expGrowthRate(r,l,t_final); % adjust growth rate to match the final length of linear growth
% growthFun = @(tnow) expGrowth(tnow,l,r); 
% .......................................
% growthFunName = 'Logistic Growth';
% kappa = linearGrowth(t_final,l,r)+0.5; %1.5; % logistic growth maximum
% r = logisticGrowthRate(kappa, r, l, t_final); % adjust growth rate to match the final length of linear growth
% growthFun = @(tnow) logisticGrowth(tnow,l,kappa,r); 
% .......................................
% growthFunName = 'Exponential Shrink';
% r = - 0.05; %expGrowthRate(r,l,t_final); % adjust growth rate to match the final length of linear growth
% l = 20.0855;
% growthFun = @(tnow) expGrowth(tnow,l,r); 
% .......................................
% growthFunName = 'Quadratic Growth';
% growthFun = @(tnow) quadraticGrowth(tnow,t_final,r,l);




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
% [u,v,epsilonu,epsilonv] = LinearFiniteDiffScheme(a,m,dt,ts,nn,ne,W,Qx,element,node,d1,d2,aa,bb,gamma,growthFun);
[u,v,epsilonu,epsilonv] = GeneralLinearFiniteDiffScheme(a,m,dt,ts,nn,ne,W,Qx,element,node,d1,d2,u0,v0,reactionFun,growthFun,noiseLevel);

toc

      
% FORMING THE SOLUTION VECTOR 
close all

figure
% --- initial state
subplot(1,2,1); hold on
plot(node,u(:,1),'r-o')
plot(node,v(:,1),'b-*')
title('initial concentration')
legend({'u','v'})
% final state
subplot(1,2,2); hold on
plot(node,u(:,ts),'r-o')
plot(node,v(:,ts),'b-*')
title('final concentration')
legend({'u','v'})

tResolution = 1000;
[uNew, xNew, tNew] = domainRemap(u,node,t,tResolution, growthFun);
[vNew, ~, ~] = domainRemap(v,node,t,tResolution, growthFun);


figure('position',[1 1 776 690],'PaperPositionMode','auto')

subplot(2,1,1)
% pcolor(t,node,u); shading interp; colorbar
pcolor(tNew,xNew,uNew); shading interp; colorbar; colormap(jet)
% caxis([0 3.5])
title('u','fontsize',14)
xlabel('time')
ylabel('x')

subplot(2,1,2)
% pcolor(t,node,v); shading interp; colorbar
pcolor(tNew,xNew,vNew); shading interp; colorbar; colormap(jet)
% caxis([0 3.5])
title('v','fontsize',14)
xlabel('time')
ylabel('x')

% ---- graphic annotation
if strcmp(growthFunName,'Logistic Growth')
    annotKappa = [', kappa = ' num2str(kappa)];
else
    annotKappa = '';
end
annotation(gcf,'textbox',...
    [0.55 0.489533011272141 0.45 0.0740740740740737],...
    'String',{[kineticFunName ' + ' growthFunName],...
              ['d_1 = ', num2str(d1), ' 	d_2 = ',num2str(d2),' 	r = ', num2str(r), annotKappa],...
              ['a = ', num2str(aa), ' 	b = ', num2str(bb), ' 	\gamma = ', num2str(gamma)] },...
    'LineStyle','none',...
    'FontWeight','bold',...
    'FontSize',10,...
    'FontName','Helvetica-Narrow',...
    'FitBoxToText','off');

% ---- save graphic
filename =['result_' datestr(now, 'yyyy-mm-dd_HH-MM-SS')];
print(['figures/' filename],'-dpng','-r0')

save(['data/', growthFunName '_' filename],'tNew','xNew','uNew','vNew')

% end
%% --- save data
save(['data/', growthFunName '_' filename],'u','v','tNew','xNew','uNew','vNew')
%% ===== plot the evolution of spatial frequency ===== %%
freq=spatialFreq(uNew,xNew(end,:));
% --- Plot 
figure('position', [464 611 747 283]);  plot(tNew,freq,'linewidth',2); %plot(uNew(:,end)); hold on;
hold on; plot(xlim,[1 1],'--r')
ylim([0.4 1.4])
xlabel('time')
ylabel('wave number')
title({'evolution of wave number on growing domain','u'})
%% ===== plot range of growing modes ===== %%
if ~exist('k','var'),k=0;end
h1=plotPeriodCountGrowth(tNew, xNew(end,:),d2,gamma,aa,bb, kineticIndex, k);
freq2n=freq.*xNew(end,:);
h2=plot(tNew, freq2n,'.r', 'displayName','actual number of periods','markersize',10);
title('Evolution of the range of allowable number of periods')
xlabel('time')
ylabel('number of periods')
legend([h1,h2],{'allowable number of periods','actual number of periods'},'location','northwest', 'fontsize',12)

%% ===== more fft ===== %%
nfft=nx;
fIdx = 1:nfft/2+1;
fuNew=fft(uNew,nfft,1);
fuNew(1,:)=0;
fuNew=sqrt(abs(fuNew(fIdx,:))); 
figure
pcolor(tNew,(fIdx-1)/l,fuNew); shading interp; colormap jet;
ylim([0 5])
ylabel('spactial frequency (cycle/unit distance)')
xlabel('time')
title('sqrt(Amplitude)')
%%
tTest=1000;
figure; plot((fIdx-1)/l,fuNew(:,tTest))
ylabel('sqrt(Amplitude)')
xlabel('spactial frequency (cycle/unit distance)')
title(['at time ' num2str(tNew(tTest))])
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