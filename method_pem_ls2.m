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
% Copyright (c) 2012 Johan Dahlin [ johan.dahlin (at) liu.se ] 
% Copyright (c) 2012 Fredrik Lindsten [ lindsten (at) isy.liu.se ]
%
% Presented at 16th IFAC Symposium on System Identification, 
% Brussels, Belgium, 2012
%
%----------------------------------------------------

function th = arls(maxOrder, data)

% Initialize estimation and validation data
data.y=double(data.y);
N = length(data.y);
ye = data.y(1:ceil(N/2));
yv = data.y(ceil(N/2)+1:end);

% Build search grid and iddata-structures
ZE = iddata(ye');
ZV = iddata(yv');
NN = (1:maxOrder)';
V = arxstruc(ZE,ZV,NN);
NN = selstruc(V,0);

% Use the System Identification toolbox to identify ARX system
sysdf=iddata(data.y');
sysarx=arx(sysdf,NN);
th = zeros(maxOrder,1);
th(1:sysarx.na) = sysarx.a(2:end);

%----------------------------------------------------
% End of File
%----------------------------------------------------
