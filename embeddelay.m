function [ex] = embeddelay(xV, m, tau)


% EMBEDDELAYS reconstructs delay trajectory from a scalar time series
% using the embedding method of delays.
% INPUTS:
%  xV    : vector of the scalar time series
%  m     : the embedding dimension
%  tau   : the delay time
% OUTPUT: 
%  ex   : Toeplitz matrix with 'm' columns and with entries the lagged
%          components of the resampled 'xV' according to the input 'tau'.

%DIMITRIADIS STAVROS 10/2007

n = length(xV);

nvec = n - (m-1)*tau;   % The length of the reconstructed set
ex = zeros(nvec,m);

for i=1:m
   ex(:,m-i+1) = xV(1+(i-1)*tau:nvec+(i-1)*tau);
end


