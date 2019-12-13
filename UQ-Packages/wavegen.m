function [etaBC] = wavegen(input,TDur,dt,Hs,Tp,f0,fHighCut,on)
%% Generate initial surface elevation based on H_S and T_P for JONSWAP
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
df=1/TDur; %frequ step

f=f0:df:fHighCut;

omega=2*pi*f;

%% Time
tSpan=0:dt:TDur;

%% Jonswap Spectrum

JS=zeros(1,length(f));

for ij=1:length(f)
   
    JS(1,ij)= jonswap(f(ij),Hs,Tp);
    
end

%% Wave generation

A=zeros(1,length(omega));
B=zeros(1,length(omega));
    
for j=1:length(omega)

    A(1,j)=sqrt(JS(1,j)*df).*input(j);
    B(1,j)=sqrt(JS(1,j)*df).*input(length(omega)+j);

end

etaBC(:,1)=sum(A(1,:).*cos(2*pi.*f.*tSpan')+B(1,:).*sin(2*pi.*f.*tSpan'),2); 

%% Making file

if on==1
    
    filename = 'initial_waves';

    Folder = pwd;

    FileName   = fullfile(Folder, sprintf('%s.iwf', filename));
    [fid, msg] = fopen(FileName, 'w');

    if fid == -1
        error('Cannot open file: %s', msg);
    end

    fprintf(fid, '# Time series generator\n');  % [EDITED]

    fprintf(fid, ['%g\n', ...
                  '%g\n'], dt,etaBC(1:end-1));  % [EDITED]
    fclose(fid);

end

return


