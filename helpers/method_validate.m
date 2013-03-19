%----------------------------------------------------
%
% Validation of model estimates
%
%----------------------------------------------------
%
% Robust ARX models with automatic order determination 
% and Student’s t innovations
%
% Authors: Johan Dahlin, Fredrik Lindsten, 
%          Thomas B. Schön, Adrian Wills.
%
% Accepted to 16th IFAC Symposium on System Identification, 
% Brussels,Belgium, 2012
%
%----------------------------------------------------

function [RMSE, VAF, yhat] = validate(A,B,y,u,y0,u0)

% Preliminaries
na = length(A);
nb = length(B);
Phi  = buildPhi2(y', u', na, nb, y0', u0');

% Calculate the predicted measurements, root mean square error
% and the fraction of variance accounted for by the estimator
% compare the R^2-measure in regression

e = y' - Phi*[A ; B];
yhat=Phi*[A ; B];
RMSE=sqrt(mean(e.^2));
VAF = 1-norm(e)/norm(y' - mean(y'));

%----------------------------------------------------
% End of File
%----------------------------------------------------
