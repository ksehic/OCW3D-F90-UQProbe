function y = runOC(BCdataX,dt)
%% Running OceanWave3D
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
%% create folder
mkdir run_case

%% Generate initial wave files
ini_wave(BCdataX,dt);

%% copy OceanWave3D file in run_case folder
copyfile ReadKinematics.m run_case
copyfile BuildStencilEven.m run_case
copyfile initial_waves.iwf run_case
copyfile beach_137 run_case
copyfile OceanWave3D.inp run_case

%% Change directory
cd run_case

%% Run OceanWave3D
command = 'OceanWave3D_devel OceanWave3D.inp';

[status,cmdout] = system(command);

%% Collect data at wind turbine
eta_ref = ReadKinematics;

%% Pick max wave height at node 100th
eta_max = max(eta_ref);

maxdata = eta_max;
y = maxdata;
%y = maxdata(:,xloc);
%% Return to original folder
cd ..

%% Delete run_case
rmdir run_case s

return