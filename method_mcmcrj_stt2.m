%----------------------------------------------------
%
% AR identification with MCMC-RJ and Student t-noise
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

function [mA N outputs]= mcmcrj_arstt(maxOrder,data)

%% Initialise
T = length(data.y);
prior.a = 0.01;
prior.b = 0.01;
proposal.ell = 0.5; % The only proposal used is for the model order jump
proposal.s2  = 0.05; % Variance for DOF random walk

% Allocate memory
miter = 4000;
A = zeros(maxOrder, miter);   % A coefficients
N = zeros(1, miter);            % Model order
S = zeros(3, miter);            % s_a^{-2}, lambda, nu

N(1) = 15; % Initialise model order
S(:,1) = [0.1 ; 0.1 ; 0.5]; % Initialise

% Build the maximal regression matrix
tmp  = buildY(data.y', maxOrder);
PHI0A = tmp(2:end,1:maxOrder);
clear tmp;

% Initialise the ARX coeffs using LS
r = N(1);
Phi = PHI0A(r:end, 1:r);
th = Phi\data.y(r+1:end)';
A(1:r,1) = th(1:r);

% We also have the latent states, but these are not logged
z = 10*ones(T-r,1); % z_{r+1:T}

%% Run loop
acceptprob = zeros(2,miter); % Log
acceptprob(:,1) = 1;

CC = 1;
for(i = 2:miter)
    if(i >= 100*CC)
        fprintf('%i : \n',i);
        CC = CC + 1;
    end
    % ---------------------- RJ-MH STEP: Sample {n,theta} ----------------------
    nap = N(i-1) + dexprnd(proposal.ell);
    % Compute acceptance probability
    if(nap > 0 && nap <= maxOrder)
        r = max(N(i-1), nap); % We remove r observations for both proposed and old orders
        z = z(r-N(i-1)+1:end); % z_{na(i-1)+1:T} --> z_{r+1:T}
        Yk = PHI0A(r:end, 1:N(1,i-1));
        Ykp = PHI0A(r:end, 1:nap);
        y = data.y(r+1:end)';
                
        sa2i = S(1,i-1); % Precision for coefficients
        Si = S(2,i-1)*z; % (Diagonal of) Precision matrix for y_{1:T}
                
        Cik = Yk'*(Yk.*repmat(Si,[1 N(i-1)])) + sa2i*eye(N(i-1));
        muk = (Cik\Yk')*(y.*Si);
        
        Cikp = Ykp'*(Ykp.*repmat(Si,[1 nap])) + sa2i*eye(nap);
        mukp = (Cikp\Ykp')*(y.*Si);
        
        logprob = 0.5*(mukp'*Cikp*mukp - muk'*Cik*muk) + 0.5*log(sa2i^nap / sa2i^N(i-1)) + 0.5*log(det(Cik)/det(Cikp));
        prob = exp(logprob);        
        acceptprob(1,i) = min(1,prob);
        
        U = rand(1);
        accept = (U < prob);
    else
        accept = false;
    end
    
    if(accept)
        N(i) = nap;
        Cikp = (Cikp+Cikp')/2;        
        tmp = mvnrnd(mukp', inv(Cikp),1)';
        A(1:nap,i) = tmp;
    else
        A(:,i) = A(:,i-1);
        N(i) = N(i-1);
    end
    % Recompute Y with new model orders (might be different than both Yk and Ykp)
    Y = PHI0A(N(1,i):end, 1:N(1,i));
    
    % GIBBS STEPS ---------------------------------------------------------
    e = data.y(N(1,i)+1:end)' - Y*A(1:N(1,i),i); % Prediction error - used below, e_{r+1:T}
    
    % Sample {z} ----------------------------------------------------------
    a = S(3,i-1)/2 + 1/2;
    b = S(3,i-1)/2 + S(2,i-1)/2*e.^2;
    % Sample from gamma distribution. N.B. inverse of b is used as parameter (scale parameter)     
    z = gamrnd(a,1./b); % z_{r+1:T}
    
    % Sample {nu, lambda, sa2i} -------------------------------------------
    % DOF, using random walk
    nu = S(3,i-1);
    nup = -1;
    while(nup < 0)
        nup = nu + sqrt(proposal.s2)*randn(1);
    end
    Ti = T-N(1,i); % Length of data record used with the current orders
    lp = -Ti*gamma(nu/2) + Ti*nu/2*log(nu/2) + nu/2*sum(log(z)) - nu/2*sum(z) + (prior.a - 1)*log(nu) - prior.b*nu;
    lpp = -Ti*gamma(nup/2) + Ti*nup/2*log(nup/2) + nup/2*sum(log(z)) - nup/2*sum(z) + (prior.a - 1)*log(nup) - prior.b*nup;
    prob = exp(lpp - lp);
    acceptprob(2,i) = min(1,prob);
    U = rand(1);
    accept = (U < prob);    
    if(accept)
        S(3,i) = nup;
    else
        S(3,i) = nu;
    end
%         Cik = Yk'*(Yk.*repmat(Si,[1 N(i-1)])) + sa2i*eye(N(i-1));
%         muk = (Cik\Yk')*(y.*Si);
%         
%         Cikp = Ykp'*(Ykp.*repmat(Si,[1 nap])) + sa2i*eye(nap);
%         mukp = (Cikp\Ykp')*(y.*Si);
%         
%         logprob = 0.5*(mukp'*Cikp*mukp - muk'*Cik*muk) + 0.5*log(sa2i^nap / sa2i^N(i-1)) + 0.5*log(det(Cik)/det(Cikp));
%         prob = exp(logprob);

    
    % Innovation precision (lambda)
    a = prior.a + 1/2*Ti;
    b = prior.b + 1/2*e'*(e.*z);
    S(2,i) = gamrnd(a,1/b);    
    
    % AR coefficient precision (delta)
    thk = A(1:N(1,i),i);;
    a = prior.a + 1/2*N(i);
    b = prior.b + 1/2*thk'*thk;
    S(1,i) = gamrnd(a,1/b);
end

%% Average
burnin=floor(miter/2);
mA=mean(A(:,burnin+1:end),2);

outputs.InPres=S(2,:);
outputs.InDOF=S(3,:);
outputs.A=A(:,burnin+1:end);
end