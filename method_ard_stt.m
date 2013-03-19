%----------------------------------------------------
%
% ARX identification with ARD priors and Student t-noise
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

function [mA mB] = method_ard_stt(maxOrder,data,robust)

%% Initialise
y = data.y(maxOrder+1:end)';
T = length(y);
prior.a = 0.01;
prior.b = 0.01;
proposal.s2  = 0.05; % Variance for DOF random walk - alternative to Stirling

% Allocate memory
miter = 4000;
A = zeros(maxOrder, miter);     % A coefficients
B = zeros(maxOrder, miter);     % B coefficients
delta = zeros(2*maxOrder, miter); % A/B precisions
S = zeros(2, miter);              % lambda, nu
E = zeros(T,miter);

% Initialise
S(:,1) = [0.1 ; 10];
delta(:,1) = 0.1;

% Build the maximal regression matrix
Phi  = buildPhi2(data.y', data.u', maxOrder, maxOrder);
Phi = Phi(maxOrder+1:end,:);

% We also have the latent states, but these are not logged
z = 10*ones(T,1); % z_{r+1:T}

%% Run loop
acceptprob = zeros(2,miter); % Log
acceptprob(:,1) = 1;

CC = 1;
for(i = 2:miter)
    if(i >= 1000*CC)
        %fprintf('%i :',i);
        CC = CC + 1;
    end

    % ---------------------- GIBBS STEP: Sample {theta} ----------------------
    Si = S(2,i-1)*z; % (Diagonal of) Precision matrix for y_{1:T}
    
    if(~robust)
        Ci = Phi'*(Phi.*repmat(Si,[1 2*maxOrder])) + diag(delta(:,i-1));
        mu = (Ci\Phi')*(y.*Si);        
        theta = mvnrnd(mu', inv((Ci+Ci')/2),1)';
    else
        [mu,R] = meancovhelper(Phi, y, Si, delta(:,i-1));
        theta = mu + R*randn(2*maxOrder,1);
    end
    A(:,i) = theta(1:maxOrder);
    B(:,i) = theta(maxOrder+1:end);

    % Compute prediction error
    e = y - Phi*theta; % Prediction error - used below
    E(:,i) = e;
    
    % ---------------------- GIBBS STEP: Sample {z} ----------------------
    a = S(2,i-1)/2 + 1/2;
    b = S(2,i-1)/2 + S(1,i-1)/2*e.^2;
    z = gamrnd(a,1./b); % Sample from gamma distribution. N.B. inverse of b is used as parameter (scale parameter)     
    
    % ---------------------- MH-STEP: Sample {nu} ----------------------
    nu = S(2,i-1);
    nup = -1;
    while(nup < 0)
        nup = nu + sqrt(proposal.s2)*randn(1);
    end    
    nup = nu + sqrt(proposal.s2)*randn(1);
    lp = -T*gamma(nu/2) + T*nu/2*log(nu/2) + nu/2*sum(log(z)) - nu/2*sum(z) + (prior.a - 1)*log(nu) - prior.b*nu;
    lpp = -T*gamma(nup/2) + T*nup/2*log(nup/2) + nup/2*sum(log(z)) - nup/2*sum(z) + (prior.a - 1)*log(nup) - prior.b*nup;
    prob = exp(lpp - lp);
    acceptprob(i) = min(1,prob);
    U = rand(1);
    accept = (U < prob);    
    if(accept)
        S(2,i) = nup;
    else
        S(2,i) = nu;
    end

    % ---------------------- GIBBS STEP: Sample {lambda,sa2i} ----------------------
    % Innovation precision
    a = prior.a + 1/2*T;
    b = prior.b + 1/2*e'*(e.*z);
    S(1,i) = gamrnd(a,1/b);
    
    % AR coefficient precision
    a = prior.a + 1/2;
    b = prior.b + 1/2*theta.^2;
    delta(:,i) = gamrnd(a,1./b);
end
      
%% Average
burnin=floor(miter/2);
mA=mean(A(:,burnin+1:end),2); mB=mean(B(:,burnin+1:end),2);
mna=maxOrder; mnb=maxOrder;
end

%----------------------------------------------------
% Help function for covariance calculation
%----------------------------------------------------
function [mu, x22, v] = meancovhelper(Phi, y, Si, delta)
[T,nth] = size(Phi);
H = Phi./repmat(delta',T,1);

M = [diag(1./Si) + Phi*H'   H ;
     H'                     diag(1./delta)];

accepted = false;
while(~accepted)
    try
        tmp = chol(M);
    catch
        warning('Regularising M');
        M = M-min(eig(M))*eye(length(M));
        tmp = chol(M);
    end
    accepted = true;
end
x11 = tmp(1:T, 1:T);
x22 = tmp(T+1:end, T+1:end);

% Compute the mean by backsubstitution
ss = x11'\y;
rr = x11\ss;
mu = Phi'*rr;
mu = mu./delta;

% Compute vector used in acceptance probability computation
v = x22'\mu;
end
%----------------------------------------------------
% End of File
%----------------------------------------------------
