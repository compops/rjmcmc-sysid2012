%----------------------------------------------------
%
% Creates estimation, validation and test data sets
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

function [data1 data2 data12 data3] = method_createdatasets(data)

% Split the data set into three parts
parts=floor(length(data.y')/3);

% Model order identification data set
data1.y=data.y(1:parts); data1.u=data.u(1:parts);

% Prediction error calculation data set
data12.y=data.y(parts+1:2*parts); data2.u=data.u(parts+1:2*parts);

% MCMC data set
data12.y=data.y(1:2*parts); data12.u=data.u(1:2*parts);

% Validation data set
data3.y=data.y(2*parts+1:end); data3.u=data.u(2*parts+1:end);

end

%----------------------------------------------------
% End of File
%----------------------------------------------------
