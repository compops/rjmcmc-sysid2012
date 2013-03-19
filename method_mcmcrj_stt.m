%----------------------------------------------------
%
% ARX identification with MCMC-RJ and Student t-noise
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

function [mA mB] = method_mcmcrj_stt(maxOrder,data,robust)

%% Initialise
T = length(data.y);
prior.a = 0.01;
prior.b = 0.01;
proposal.ell = 0.5; % The only proposal used is for the model order jump
proposal.s2  = 0.05; % Variance for DOF random walk

% Allocate memory
miter = 4000;
A = zeros(maxOrder, miter);   % A coefficients
B = zeros(maxOrder, miter);   % B coefficients
N = zeros(2, miter);            % Model order
S = zeros(3, miter);            % s_a^{-2}, lambda, nu

N(:,1) = 15; % Initialise model order
S(:,1) = [0.1 ; 0.1 ; 0.5]; % Initialise

% Build the maximal regression matrix
tmp  = buildPhi2(data.y', data.u', maxOrder, maxOrder);
PHI0A = tmp(2:end,1:maxOrder);
PHI0B = tmp(2:end,maxOrder+1:end);
clear tmp;

% Initialise the ARX coeffs using LS
r = N(1,1);
s = N(2,1);
Phi = [PHI0A(r:end, 1:r), PHI0B(r:end, 1:s)];
th = Phi\data.y(r+1:end)';
A(1:r,1) = th(1:r);
B(1:s,1) = th(r+1:end);

% We also have the latent states, but these are not logged
z = 10*ones(T-r,1); % z_{r+1:T}

%% Run loop
acceptprob = zeros(2,miter); % Log
acceptprob(:,1) = 1;

CC = 1;
for(i = 2:miter)
    if(i >= 1000*CC)
        %fprintf('%i :',i);
        CC = CC + 1;
    end
    % ---------------------- RJ-MH STEP: Sample {n,theta} ----------------------
    nap = N(1,i-1) + dexprnd(proposal.ell);
    nbp = N(2,i-1) + dexprnd(proposal.ell);

    % Compute acceptance probability
    if(nap > 0 && nap <= maxOrder && nbp > 0 && nbp <= maxOrder);
        r = max(N(1,i-1), nap); % We remove r observations for both proposed and old orders
        z = z(r-N(1,i-1)+1:end); % z_{na(i-1)+1:T} --> z_{r+1:T}
        Yk = [PHI0A(r:end, 1:N(1,i-1)) PHI0B(r:end, 1:N(2,i-1))];
        Ykp = [PHI0A(r:end, 1:nap) PHI0B(r:end, 1:nbp)];
        y = data.y(r+1:end)';
                
        sa2i = S(1,i-1); % Precision for coefficients
        Si = S(2,i-1)*z; % (Diagonal of) Precision matrix for y_{1:T}
        
        if(~robust)        
            Cik = Yk'*(Yk.*repmat(Si,[1 (N(1,i-1)+N(2,i-1))])) + sa2i*eye(N(1,i-1)+N(2,i-1));
            muk = (Cik\Yk')*(y.*Si);
            
            Cikp = Ykp'*(Ykp.*repmat(Si,[1 nap+nbp])) + sa2i*eye(nap+nbp);
            mukp = (Cikp\Ykp')*(y.*Si);
            
            logprob = 0.5*(mukp'*Cikp*mukp - muk'*Cik*muk) + 0.5*log(sa2i^(nap+nbp) / sa2i^(N(1,i-1)+N(2,i-1))) + 0.5*log(det(Cik)/det(Cikp));
            prob = exp(logprob);
        else
            [muk,Rk,vk] = meancovhelper(Yk, y, Si, sa2i);
            [mukp,Rkp,vkp] = meancovhelper(Ykp, y, Si, sa2i);
            
            logprob = 0.5*(vkp'*vkp - vk'*vk) + 0.5*(log(sa2i^(nap+nbp)) - log(sa2i^(N(1,i-1)+N(2,i-1)))) + (log(det(Rkp))-log(det(Rk)));
            prob = exp(logprob);
        end
        
        acceptprob(1,i) = min(1,prob);
        
        U = rand(1);
        accept = (U < prob);
    else
        accept = false;
    end
    
    if(accept)
        N(:,i) = [nap ; nbp];
        if(~robust)
            Cikp = (Cikp+Cikp')/2;        
            tmp = mvnrnd(mukp', inv(Cikp),1)';
        else
            tmp = mukp + Rkp*randn(nap+nbp,1);
        end
        A(1:nap,i) = tmp(1:nap);
        B(1:nbp,i) = tmp(nap+1:end);
    else
        A(:,i) = A(:,i-1);
        B(:,i) = B(:,i-1);
        N(:,i) = N(:,i-1);
    end

    % Recompute Y with new model orders (might be different than both Yk and Ykp)
    Y = [PHI0A(N(1,i):end, 1:N(1,i)) PHI0B(N(1,i):end, 1:N(2,i))];
   
    % Prediction error - used below, e_{r+1:T}
    e = data.y(N(1,i)+1:end)' - Y*[A(1:N(1,i),i) ; B(1:N(2,i),i)]; 
    
    % ---------------------- GIBBS STEP: Sample {z} ----------------------
    a = S(3,i-1)/2 + 1/2;
    b = S(3,i-1)/2 + S(2,i-1)/2*e.^2;

    % Sample from gamma distribution. N.B. inverse of b is used as parameter (scale parameter)     
    z = gamrnd(a,1./b); % z_{r+1:T}

    % ---------------------- MH-STEP: Sample {z} ----------------------
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

    % ---------------------- GIBBS STEP: Sample {lamda,sa2i} ----------------------
    % Innovation precision (lambda)
    a = prior.a + 1/2*Ti;
    b = prior.b + 1/2*e'*(e.*z);
    S(2,i) = gamrnd(a,1/b);    
    
    % AR coefficient precision (delta)
    thk = [A(1:N(1,i),i); B(1:N(2,i),i)];
    a = prior.a + 1/2*(N(1,i) + N(2,i));
    b = prior.b + 1/2*thk'*thk;
    S(1,i) = gamrnd(a,1/b);
end

%% Average
burnin=floor(miter/2);
mA=mean(A(:,burnin+1:end),2); mB=mean(B(:,burnin+1:end),2);

end

%----------------------------------------------------
% Help function for covariance calculation
%----------------------------------------------------
function [mu, x22, v] = meancovhelper(Phi, y, Si, delta)
[T,nth] = size(Phi);

M = [diag(1./Si) + Phi*Phi'/delta   Phi/delta ;
     Phi'/delta                     eye(nth)/delta];

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
mu = Phi'*rr/delta;

% Compute vector used in acceptance probability computation
v = x22'\mu;

end
%----------------------------------------------------
% End of File
%----------------------------------------------------
