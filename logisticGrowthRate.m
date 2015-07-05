function logisticRate = logisticGrowthRate( kappa, linearRate, L0, t_final )
%LOGISTICGROWTHRATE recalculate the growth rate for logistic growth so that
%the final length of the domain is the same with linear growth for the same
%growing time.
%   Detailed explanation goes here

logisticRate = 1/t_final*log((L0+linearRate*t_final)/(L0-linearRate*t_final/(kappa-1)));
end

