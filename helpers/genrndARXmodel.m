%----------------------------------------------------
%
% Generation of random ARX systems with 
% outliers/different noise distributions
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

function [sys0, data] = genrndARXmodel(na,nb,T,distribution)
% genrndARXmodel Simulates a random ARX-system
% [SYS0, DATA, DATA1, DATA2, DATA12, DATA3] = 
%                   genrndARXmodel(NA,NB,T,distribution) 
% returns the system struct (SYS0) containing the model order and the model 
% coefficents, and simulated random data (DATA) u, y and e for the system. 
% The data set in splitted into three (equally large) parts: DATA1, DATA2,
% and DATA3. DATA12 contain the first to parts together. (Training data)
%
% DISTRIBUTION determines the pdf from which the innovations are generated
% gaussian or sst (student-t). If 'deterministic' noise is specified then 
% random outliers are added to the output signal. If 'sstnoise' is 
% specified then student-t distributed noise is added to the output signal.
%
% If only NA is provided the model orders are randomly selected with NA as 
% the maximum possible order. Providing NA and NB yields a system of those 
% model orders. T is the number of time steps that are to be simulated
%

% Check number of inputs are rearrange of use defaults.
if (nargin < 4) 
   if (nargin == 2); distribution='gaussian'; nmax=na; T=nb; end
   if (nargin == 3); nmax=na; distribution=T; T=nb; end
   
   nb=2; na=1; 
   
    % Check for proper system
    while ((na < nb))
        na = randi(nmax);
        nb = randi(nmax);
    end
end

% Generate more data than needed to obtain stationary process.
T=2*T;

%% Generate noise model

% Generate the nominator polynomial
nom = genpolynomial(nb);

% Generate the denominator polynomial (requiring proper TF)
denom = genpolynomial(na);

sys0.n = [na nb]; sys0.type = 'arx stt';

a = poly(denom.roots); b = poly(nom.roots);
sys0.a = a(2:end); sys0.b = b(1:nb);

theta = [sys0.a sys0.b]';

%% Input
u = [zeros(1,nb+T/2) 0.1*randn(1,T/2)]; % u_t = u(t+nb)

%% Noise
sigma = 1;        % Gaussian noise variance
lambda = 1;       % Innovation precision
nu = 2;           % Innovation degrees of freedom
outliersize = 5;  % Size of outliers
lambdanoise = 10;  % Innovation precision of sttnoise
nunoise = 2;      % Innovation degrees of freedom of stt noise

% Initialize
data.tVec = (1:T/2);
y = zeros(1, T+na);

% Generate Student-t Noise
if strcmp(distribution,'stt'),
    e = 1/sqrt(lambda)*trnd(nu, 1, T);
    sys0.lambda=lambda; sys0.nu=nu;
end

% Generate Gaussian Noise
if ~strcmp(distribution,'stt'),
    e = 1/sqrt(lambda)*randn(1, T);
    sys0.noisevariance=1/lambda;
end

for(t = 1:T)
    phi = [u(t:t+nb-1) -y(t:t+na-1)]; % u(t-nb) .. u(t-1), -y(t-na) .. -y(t-1)
    phi = phi((na+nb):-1:1); % Flip
    y(t+na) = phi*theta + e(t);
end

data.y = y(na+1+T/2:end); % discard the na initial zeros
data.u = u(nb+1+T/2:end); % discard the nb initial zeros
data.e = e(1+T/2:end);
T = T/2;

%% Add outliers to output signal
if strcmp(distribution,'sttnoise'),
    outputnoise = 1/sqrt(lambdanoise)*trnd(nunoise, 1, T);    
    % Add outliers    
    data.y=data.y+outputnoise;
end

if strcmp(distribution,'deterministic'),
    Te = floor(2*T/3)-30;
    %noutliers=randi(floor(0.01*Te));
    %toutliers=randsample(1:length(e),randi(length(e)),0);
    noutliers = round((1+rand(1)*2)*Te/100); % 1-3 % outliers
    toutliers=randsample(Te,noutliers,0);
    
    % Store information
    sys0.noutliers=noutliers; 
    sys0.sizeoutliers=outliersize*max(abs(y));
    sys0.toutliers = toutliers;
    
    % Add outliers    
    data.y(toutliers)=data.y(toutliers)+sys0.sizeoutliers*(2*rand(1,noutliers)-1);
elseif strcmp(distribution,'missingdata'),
    missingnoise = 0.1; % Standard deviation of the noise used when we miss observation
    Te = floor(2*T/3)-30;
    %Te = T;
    %noutliers=max(1,randi(floor(0.05*Te)));
%     noutliers = round(0.02*Te); % 2 % outliers
    noutliers = round((1+rand(1)*2)*Te/100); % 1-3 % outliers
    toutliers=randsample(Te,noutliers,0);
    
    % Store information
    sys0.noutliers=noutliers;     
    sys0.toutliers = toutliers;
    
    % Add outliers    
    data.y(toutliers) = missingnoise*randn(1,noutliers);
end


%----------------------------------------------------
% Generate random polynomial with stable roots
%----------------------------------------------------
function polyn = genpolynomial(N)
% genpolynomial Generate a random polynomial with N roots.
% [polyn] = genpolynomial(N) returns a struct with the model order (n),
% the number of complex roots (k) and real roots (l), the roots are also
% returned in the struct (roots). The number of complex and real roots 
% are randomly selected.
%

polyn.n=N; % number of roots
polyn.k=2*floor(randi(polyn.n)/2); % number of complex roots
polyn.l=polyn.n-polyn.k; % number of real roots

% Generate real roots
polyn.roots=-0.95+rand(polyn.l,1)*(0.95+0.95);

% Generate complex roots
for (j=1:polyn.k/2)
    
    correct=0;
    while (~correct)    
       % Generate point in upper unit square
       x1=-0.95+rand*(0.95+0.95); 
       x2=rand*0.95; 
    
        %Accept if within the unit circle
        if (sqrt(x1.^2+x2.^2) < 0.95)
            correct=1;
        end
    end
    
    % Store root with its complex conjugate
    polyn.roots=[polyn.roots; x1+x2*sqrt(-1); x1-x2*sqrt(-1)];
end
%----------------------------------------------------
% End of File
%----------------------------------------------------
