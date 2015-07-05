function expRate = expGrowthRate( linearRate, L0, t_final )
%EXPGROWTHRATE recalculate the growth rate for expoential growth so that
%the final length of the domain is the same with linear growth for the same
%growing time.
%   usage: 

expRate = log(1+linearRate*t_final/L0)/t_final;
end

