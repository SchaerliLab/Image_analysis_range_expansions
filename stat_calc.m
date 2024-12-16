function [weightedmean,confintcalced,weights] = stat_calc(a,w)
% a: the data poins
% w: std of the data points
% weightedmean: weighted mean with Bisquare weights
% confintcalced: calculated 95 % confidence interval
% weights: bisquare weights with this method: https://uk.mathworks.com/help/stats/robustfit.html#mw_84ae84c6-63ed-4c22-93bb-4d0064b07c51

if any(isnan(a),'all')
    weightedmean = nan;
    confintcalced = nan;
    warning('nan value is stat_calc')
    return
else

% fitting a constant line with weights and calulating the confidence
% intervals
fun = @(k,x) k+0*x;

fit_options = fitoptions(fun,'StartPoint',mean(a,'all'),'Robust','Bisquare');
if nargin==2
    fit_options.Weights = 1./w;
end

fitobject = fit((1:numel(a)).',a(:),fun,fit_options);
ci = confint(fitobject);

weightedmean = fitobject.k;
confintcalced = fitobject.k-ci(1);

% calculate the bisquare weights for the colonies
% https://uk.mathworks.com/help/stats/robustfit.html#mw_84ae84c6-63ed-4c22-93bb-4d0064b07c51
% get the residuals
% residuals = output.residuals;
residuals = a-fitobject.k;
% the vector of leverage values
h = 0;%leverage(x);
% median absolute deviation
MAD_calc = mad(residuals, 1);
% an estimate of the standard deviation of the error term
ss = MAD_calc/0.6745; 
% tuning constant
tune = 4.685;

r = residuals ./ (tune * ss * sqrt(1 - h));
% the calculated bisquared weights
weights = (abs(r) < 1) .* (1 - r.^2).^2;

end
end

