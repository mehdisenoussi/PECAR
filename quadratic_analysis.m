function [p1,p2,delta] = quadratic_analysis(pBoth,pNone)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This function calculates for each individual p1 and p2 (i.e. the amount %%%
%%% of attention at each of the two probe locations) based on pBoth         %%%
%%% (probability to report correctly the two probes) and pNone (probability %%%
%%% to report correctly none of the probes).                                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INPUT:
% pBOTH and pNONE must be vectors of the averaged probabilities for each probe delay
%%% OUTPUT:
% p1 and p2 must be vectors of p1 and p2 values for each probe delay

%% Computation of the quadratic parameters assuming p1 and p2 are independent
% Basic quadratic theory says: x1+x2=sigma and x1*x2=phi when X^2 - sigma*X + phi = 0
% After simplification, we can say:

sigma = 1+pBoth-pNone;
phi = pBoth;

%% Computation of the quadratic discriminant (delta)
delta = sigma.^2 - 4.*phi;

%% Keep imaginary delta values (DO NOT reset to zero or flip imaginary into positive values)
holder = delta<0;
delta_sq = delta;
delta_sq(holder) = -delta_sq(holder);
delta_sq = sqrt(delta_sq);
delta_sq(holder) = -delta_sq(holder);

%% Compute P1 and P2: P1 is either higher or equal to P2 by mathematical construction
% P1 and P2 are the two solutions of the following second-degree equation:
% X^2 - sigma*X + phi = 0
p1 = (sigma + delta_sq)./2;
p2 = (sigma - delta_sq)./2;

end