%----------------------------------------------------
%
% Identification of ARX process using PEM (Least Squares)
%
%----------------------------------------------------
%
% Hierarchical Bayesian ARX models for robust inference
%
% Authors: Johan Dahlin, Fredrik Lindsten, 
%          Thomas B. Schön, Adrian Wills.
%
% Presented at 16th IFAC Symposium on System Identification, 
% Brussels, Belgium, 2012
%
%----------------------------------------------------

function [A B] = method_pem_ls(maxOrder, data)
% [A B] = METHOD_PREM_LS(maxOrder, data)
% Inputs: maxOrder (maximum order of the search grid)
%         data (data divided into estimation and validation set
% Output: A and B polynoms.
%

% Initialize estimation and validation data
N = length(data.y);
ye = data.y(1:ceil(N/2)); ue = data.u(1:ceil(N/2));
yv = data.y(ceil(N/2)+1:end); uv = data.u(ceil(N/2)+1:end);

% Build search grid and iddata-structures
ZE = iddata(ye',ue');
ZV = iddata(yv',uv');
[X, Y] = meshgrid(1:maxOrder);
NN = [X(:) Y(:) ones(maxOrder^2,1)];
V = arxstruc(ZE,ZV,NN);
NN = selstruc(V,0);

% Use the System Identification toolbox to identify ARX system
sysdf=iddata(data.y',data.u');
sysarx=arx(sysdf,NN);
A = zeros(maxOrder,1);
A(1:sysarx.na) = sysarx.a(2:end);
B = zeros(maxOrder,1);
B(1:sysarx.nb) = sysarx.b(2:end); % 1 step delay

%----------------------------------------------------
% End of File
%----------------------------------------------------
