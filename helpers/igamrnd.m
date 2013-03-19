function R = igamrnd(a,b,varargin)
% IGAMRND Random arrays from inverse gamma distribution.
%    R = IGAMRND(A,B) returns an array of random numbers chosen from the
%    inverse gamma distribution with shape parameter A and scale parameter B.
% 
%    R = GAMRND(A,B,M,N,...) or R = GAMRND(A,B,[M,N,...]) returns an
%    M-by-N-by-... array.

R = 1./gamrnd(a,1/b,varargin{:});