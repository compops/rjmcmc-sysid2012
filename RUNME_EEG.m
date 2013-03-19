%----------------------------------------------------
%
% Example of EEG data analysis
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

%% Initialization
load eegdataforpublication
N=length(y);

%% Analysis
% RJMCMC step - inference and validation
Nval=floor(2*N/3); y12 = y(1:Nval); data12.y = y12; y3 = y(Nval+1:end);
[th tmp outputs]= method_mcmcrj_stt2(30,data12); Norder=tmp;
[rmse, vaf, yhat] = method_validate2(th, y3, y12);

% LS - inference and validation
thLS = method_pem_ls2(30,data12);
[rmseLS, vafLS, yhatLS] = method_validate2(thLS, y3, y12);        


%% Plotting
figure();
subplot(2,1,1); plot(1:N,y,Nval+1:N,yhat,'r',Nval+1:N,yhatLS,'g');
xlabel('time (t)'); ylabel('signal (y)');

subplot(2,2,3);
burnin=1000; miter=4000;
d = histc(Norder(burnin+1:end),(1:30)-0.5);
stem(1:30,d/(miter-burnin),'color','blue','linewidth',2)
xlabel('model order (k)'); ylabel('probability');
axis([0 20 0 0.6]);

subplot(2,2,4)
qqplot(y); title(''); ylabel('Sample quantiles');

%% Optional plot of the DOF in the Student-t distribution
subplot(2,2,4)
DOF=outputs.InDOF(burnin+1:end);
[bwd den xmesh] = kde(DOF);
plot(xmesh,den,'color','blue','linewidth',2);
title('Innovation DOF (nu)');
axis([2.5 4.5 0 6]);






