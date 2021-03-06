%returns X, Y, & R^2 Values for a linear fit
function [RSQ P] = oneOneLinearFit(X,Y)
% Remove any present NaN values
X = X(isfinite(X(:,1)),:);
Y = Y(isfinite(Y(:,1)),:);
% Create 1000 point scale for span of X and perform fit
XX = min(X):max(X)/1000:max(X);
%[P] = polyfit(X,Y,1);
YY = XX.*1;%+P(2);
P(1) = 1;
P(2) = 0;
% Evaluate correlation coefficient
Poly = polyval(P,X);
Resid = Y - Poly;
ResidSS = sum(Resid.^2);
ResidTotal = (length(Y)-1)*var(Y);
RSQ = 1 - ResidSS/ResidTotal;
end