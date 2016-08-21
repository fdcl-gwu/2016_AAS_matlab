% 4 Mar 2016
% test speed improvements of polyhedron potential code
clc
clear all
close all

constants = load_constants('castalia');

% perform preliminary computations for polyhedron potential model
num_runs = 100;
ef_time = zeros(num_runs,1);
pot_time = zeros(num_runs,1);
mex_time = zeros(num_runs,1);

rv = [-0.4358   -0.4815   -0.3051];
F = [0.1527   -0.1658   -0.3192;
   -0.1658    0.1800    0.3466;
   -0.3192    0.3466    0.6673];
wface = 7.8289e-04;


% for ii = 1:num_runs
%     start = tic;
%     asteroid_grav = polyhedron_shape_input(constants);
%     ef_time(ii) = toc(start);
%     clear asteroid_grav
% end
% fprintf('Mean E,F calculation %3.6e seconds\n',mean(ef_time));

asteroid_grav = polyhedron_shape_input(constants);
% calculate potential at a single state
state = [1;0.2;0];
for ii = 1:num_runs
    
    start = tic;
    [U,Ug,Ug_mat, Ulap] = polyhedron_potential(state, asteroid_grav);
    pot_time(ii) = toc(start);
    start = tic;
    [U,Ug,Ug_mat, Ulap] = polyhedron_potential_mex(state, asteroid_grav);
    mex_time(ii) = toc(start);
end

fprintf('Mean potential matlab %3.6e seconds\n',mean(pot_time));
fprintf('Mean potential MEX %3.6e seconds\n',mean(mex_time));

mat_time = zeros(num_runs,1);
codegen_time = zeros(num_runs,1);
eigen_time = zeros(num_runs,1);

for ii = 1:num_runs
    rv = rand(1,3);
    F = rand(3,3);
    wface = rand(1,1);
    
    start = tic;
    [U_face, U_grad_face, U_grad_mat_face] = face_contribution(rv,F,wface);
    mat_time(ii) = toc(start);
    
    start = tic;
    [U_face, U_grad_face, U_grad_mat_face] = face_contribution_mex(rv,F,wface);
    codegen_time(ii) = toc(start);
    
    start = tic;
    [U_face, U_grad_face, U_grad_mat_face] = face_contribution_eigen(rv,F,wface);
    eigen_time(ii) = toc(start);
    
end

fprintf('Matlab %3.6e seconds\n',mean(mat_time));
fprintf('Codegen MEX %3.6e seconds\n',mean(codegen_time));
fprintf('Eigen MEX %3.6e seconds\n',mean(eigen_time));