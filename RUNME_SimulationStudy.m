%----------------------------------------------------
%
% Example of Simulation experiment
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

%% Initialization
clear

% maxOrder - the maximum order of the generating polynomials (ARX order)
% T - the number of observations
% Ntrials - the number of generated systems
maxOrder=30; T=450; Ntrials=10; 
plotON = true; modeltype = 'stt';

% Preallocate outputs
LSrmse=zeros(Ntrials,1); MCMCrmse=zeros(Ntrials,1); MCMCrmse2=zeros(Ntrials,1);
LSvaf=zeros(Ntrials,1); MCMCvaf=zeros(Ntrials,1); MCMCvaf2=zeros(Ntrials,1);
mALS=zeros(maxOrder,Ntrials); mBLS=zeros(maxOrder,Ntrials);
mA=zeros(maxOrder,Ntrials); mB=zeros(maxOrder,Ntrials);
mA2=zeros(maxOrder,Ntrials); mB2=zeros(maxOrder,Ntrials);
yhatLS=zeros(ceil(T/3),Ntrials); yhatMCMCJ=zeros(ceil(T/3),Ntrials);
yhatMCMCARD=zeros(ceil(T/3),Ntrials);


% Main loop, for Ntrials systems
for i=1:Ntrials
     fprintf('%i :',i);
        
     % Generate data
     marker=0;
     while (marker==0)
         [sys0(i), data(i)] = genrndARXmodel(maxOrder,T,modeltype);
         % Check for systems with very large values (this will give us problem
         % with some matrix inversions in the covariance estimation).
         if (max(data(i).y)<1000); marker=1; end
     end
        
     % Divide the data into validation and estimation sests
     [data1(i), data2(i), data12(i), data3(i)] = method_createdatasets(data(i));
     fprintf('.');
	
     % Prediction Error Metod (Least Squares)
     [mALS(:,i) mBLS(:,i)] = method_pem_ls(maxOrder, data12(i));
     [LSrmse(i), LSvaf(i), yhatLS(:,i)] = ...
        method_validate(mALS(:,i), mBLS(:,i), data3(i).y, data3(i).u, data12(i).y, data12(i).u);
     fprintf('.');
	
     % RJ-MCMC and ARD prior methods
     [mA(:,i) mB(:,i)] = method_mcmcrj_stt(maxOrder,data12(i),0);
     [mA2(:,i) mB2(:,i)] = method_ard_stt(maxOrder,data12(i),0);
     [MCMCrmse(i), MCMCvaf(i), yhatMCMCJ(:,i)] = ...
         method_validate(mA(:,i), mB(:,i), data3(i).y, data3(i).u, data12(i).y, data12(i).u);
     [MCMCrmse2(i), MCMCvaf2(i), yhatMCMCARD(:,i)] = ...
         method_validate(mA2(:,i), mB2(:,i), data3(i).y, data3(i).u, data12(i).y, data12(i).u);
    fprintf('.\n');

end

%% Plotting
if(plotON)
	figure(22);
	subplot1(1,2,'Gap',[0.06 0.06])
	subplot1(1)
	plot(1:Ntrials,(MCMCvaf-LSvaf),'black','linewidth',1);
	hold on; plot([1 Ntrials],[0 0],'--','color','black'); hold off
    title('RJ-MCMC vs PEM','interpreter','latex','FontSize',17);
    xlabel('system','interpreter','latex','FontSize',17);
	ylabel('model fit difference','interpreter','latex','FontSize',17);
	set(gca,'FontSize',14)

    subplot1(2)
    plot(1:Ntrials,(MCMCvaf2-LSvaf),'black','linewidth',1);
	hold on; plot([1 Ntrials],[0 0],'--','color','black'); hold off
	title('ARD prior vs PEM','interpreter','latex','FontSize',17);
	xlabel('system','interpreter','latex','FontSize',17);
	ylabel('model fit difference','interpreter','latex','FontSize',17);
	set(gca,'FontSize',14)

end


end





