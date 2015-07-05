aa = 0.1;     bb = 1;


periods = 1:10; n_periods= length(periods);
diffusionRatio = 10:10:100; n_diffustionRation = length(diffusionRatio);

relation = zeros(n_diffustionRation, n_periods);
for n=1:n_periods
    for d = 1:n_diffustionRation
        relation(d,n) = gammaByHalfWaves(2*periods(n),aa,bb,diffusionRatio(d),1,2);
    end
end

figure; hold on
cl=hsv(n_periods);
for i=1:n_periods
    ph(i)=plot(diffusionRatio(:), relation(:,i),'color',cl(i,:)); 
end
set(gca,'color','k','linewidth',2)
title('relationship between d and \gamma for fixed fastest growing wave mode')
xlabel('d')
ylabel('\gamma')
legend([ph(1),ph(end)],...
    {['spatial frequency = ', num2str(periods(1))],['spatial frequency = ' num2str(periods(n_periods))]},...
    'textcolor','w','fontsize',14)