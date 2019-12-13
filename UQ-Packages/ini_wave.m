function [] = ini_wave(etaBC,dt)
%% Generate initial wave file
%{
---------------------------------------------------------------------------
Created by:
Kenan Šehić (kense@dtu.dk; kenosehic@gmail.com)
Department of Applied Mathematics and Computer Science
Technical University of Denmark
Licence: Copyright (C) 2019 Kenan Šehić DTU Compute, Technical University of Denmark

Cite: Šehić K., Bredmose H., Sørensen J.D., Karamehmedović M.: Low-dimensional representation of wave generation to quantify extreme events, TBD
---------------------------------------------------------------------------
Version December 2019
---------------------------------------------------------------------------
%}
%% Procedure

fprintf('Making Initial Waves...\n');

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

fprintf('Done!\n');

return
