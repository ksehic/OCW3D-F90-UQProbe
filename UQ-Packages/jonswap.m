function s=jonswap(f,Hs,Tp)
%% JONSWAP spectrum
%{
---------------------------------------------------------------------------
Created by:
Kenan Šehić (kense@dtu.dk; kenosehic@gmail.com)
Department of Applied Mathematics and Computer Science
Technical University of Denmark
Licence: Copyright (C) 2019 Kenan Šehić DTU Compute, Technical University of Denmark

Cite: Šehić K., Bredmose H., Sørensen J.D., Karamehmedović M.: Low-dimensional representation of wave generation to quantify extreme events, TBD
Status: Submitted - Journal of Engineering Mathematics Dec 2019
---------------------------------------------------------------------------
Version December 2019
---------------------------------------------------------------------------
%}
%% Procedure
fp=1/Tp; % Peak frequency

TH=Tp/sqrt(Hs);

%peak-shape parameter
if TH<=3.6
    gamma=5;
elseif (3.6<TH) && (TH<=5)
    gamma=exp(5.75-1.15*TH);
else
    gamma=1;
end

if f<=fp
    sigma=0.07;
else
    sigma=0.09;
end

s=0.3125*Hs^2*Tp*(f/fp)^(-5)*exp(-1.25*(f/fp)^(-4))*(1-0.287*log(gamma))*gamma^(exp(-0.5*((f/fp-1)/sigma)^2));

return
