function [L,w] = weighted_aic(delta)

%  INPUT:  delta = vector of AIC differences
%
%  OUTPUT:     L = vector of relative likelihood estimates
%              w = vector of weighted AIC values
% 
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %
%                                                                         %
%  Computes weighted AIC value for model selection                        %
%  See Wagenmakers & Farrel (2004)                                        %
%                                                                         %
%  Last updated: 04.14.2020                                               %
%  Jonathan K Doyon, PhD                                                  %
%  jdoyon@gwu.edu                                                         %
%                                                                         %
%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %


% estimate relative likelihood L of model i
for i = 1:length(delta)
    L(i,1) = exp( -0.5 * delta(i) );
end
% normalize likelihoods to obtain weighted AIC
for i = 1:length(delta)
    w(i,1) = L(i) / sum(L);
end
end

