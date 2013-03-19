function [RMSE, VAF, yhat] = validate_ar(th,y,y0)
% y = Phi*th

na = length(th);
Phi  = buildY([y0' ; y'],na);
Phi = Phi(length(y0)+1:end,:);
yhat=Phi*th;
e = y' - yhat;
RMSE=sqrt(mean(e.^2));
%VAF=1-var(e)/var(y');
err = norm(e);
meanerr = norm(y' - mean(y'));
VAF = 1-err/meanerr;
