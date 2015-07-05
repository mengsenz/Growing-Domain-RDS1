function hp=plotPeriodCountGrowth( t, Lt, d, gamma,a,b, kineticModel, kGM)
%PLOTPERIODCOUNTGROWTH Summary of this function goes here
%   usage:  plotPeriodCountGrowth( t, Lt, d, gamma,a,b, kineticModel, kGM )
if nargin<8, kGM=0; end % Gierer-Meinhardt saturation
if nargin<7, kineticModel = 1; end 

% nt=length(Lt);
initRange = periodCountRange(1,d,gamma,a,b, kineticModel, kGM);
n_max = initRange(1)*Lt(:);
n_min = initRange(2)*Lt(:);
upperbound=max(n_max);


figure('Position',[485 256 736 675])
h=area(t,[n_min,n_max-n_min]); hold on
try
    h(1).FaceColor='w';
    h(2).FaceColor=0.941*ones(1,3);
catch ME
    set(h(1),'FaceColor','w')
    set(h(2),'FaceColor',0.941*ones(1,3))
end

[X,Y]=meshgrid(t,0:0.5:max(n_max));
X=num2cell(X,1);
Y=num2cell(Y,1);
n_max = num2cell(n_max);
n_min = num2cell(n_min);
YInRange = cellfun(@(x,xmin,xmax) x(x>=xmin & x<=xmax), Y, n_min', n_max','uniformoutput',0); 
XInRange = cellfun(@(y,x,xmin,xmax) y(x>=xmin & x<=xmax), X, Y, n_min', n_max','uniformoutput',0); 

XInRange = cell2mat(XInRange');
YInRange = cell2mat(YInRange');

hp=plot(XInRange,YInRange, '.','color', [0.1 0.447 0.741],'markersize',0.5,'DisplayName','allowable number of periods');
ylim([0,upperbound])
end

