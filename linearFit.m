%returns X, Y, & R^2 Values for a linear fit
function [XX, YY, RSQ, PP] = linearFit(X,Y)
% Remove any present NaN values
X = X(isfinite(X(:,1)),:);
Y = Y(isfinite(Y(:,1)),:);
% Create 1000 point scale for span of X and perform fit
XX = min(X):max(X)/1000:max(X);
% Check for wrong vector direction
try
    [PP] = polyfit(X,Y,1);
catch
    X = X.';
    [PP] = polyfit(X,Y,1);
end
YY = XX.*PP(1)+PP(2);
% Evaluate correlation coefficient
Poly = polyval(PP,X);
Resid = Y - Poly;
ResidSS = sum(Resid.^2);
ResidTotal = (length(Y)-1)*var(Y);
RSQ = 1 - ResidSS/ResidTotal;
end