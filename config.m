function [fd, params, dof] = config()
% Config: Configuration file for the phase field model
%   [fd, params, dof] = config()
%
% Inputs:
%   USER DEFINED PARAMETERS
%
% Outputs:
%   fd:     Finite difference matrices and discretization
%   params: Parameters for the problem
%   dof:    Degrees of freedom for the problem
%
% By Naghmeh Akhavan March 2025

% Set the parameters for the problem

params = struct();  % Initialize the parameters struct
dof = struct();     % Initialize the degrees of freedom struct
fd = struct();      % Initialize the finite difference struct


%% Biological parameters
% Cell type I (Nurse cells)
params.cell_radi_nurse = [0.0005, 0.0005, 0.0005, ...
                            0.0005, 0.0005, 0.0005];
params.x_centers_nurse = [1.5, 2.3, 3.1, ...
                            1.4, 2.3, 3.0];
params.y_centers_nurse = [1.6, 1.9, 2., ...
                            3.3, 3.1, 3];
params.target_volumes_nurse = [0.5, 0.6, 1., ...
                                0.6, 0.8, 1.1];
% Cell type II (Oocyte)
params.cell_radi_oocyte = 0.0005;
params.x_centers_oocyte = 3.4;
params.y_centers_oocyte = 2.5;
params.target_volumes_oocyte = 1.95;
% Cell type III (Cluster)
params.cell_radi_cluster = 0.0005;
params.x_centers_cluster = .95;
params.y_centers_cluster = 2.5;
params.target_volumes_cluster = 0.2;

params.dval = [0.001; 0.001; 0.0005];
params.alpha = 100;
params.beta_s = .9;        % epithelial repulsion
params.eta_s = 0.007;       % epithelial adhesion
params.gamma = 0.01;        % gamma > eta
params.beta = [0.25, 0.25, 0.25;  %NC, Oocy, Clus
                0.25, 0, 0.3;
                0.25, 0.3, 0.];
params.eta = [0.003, 0.004, 0.008;   %NC, Oocy, Clus
                0.004, 0, 0.005;
                0.008, 0.005, 0];

%% Discretization parameters
params.x_length = 5;            % Length of the domain
params.y_length = 5;            % Width of the domain
params.h_space = 0.05;          % space mesh size
params.h_time = 0.05;          % time step size
params.t_final = 500;           % final time
params.tol = 1e-12;            % tolerance for the solver

dof.n_saved = 500;  % Number of solutions saved
dof.n_saved_c = 500;  % Number of solutions saved for concentration
dof.n_saved_y = 500;  % Number of solutions saved for concentration

%% Boolean parameters
params.run_steadyQ = 1;  % Run the problem without concentration (steady state)
params.run_concenQ = 1;  % Run the concentration problem

params.run_tensionQ =1;  % Run the problem with tension force
params.fixed_ConcenQ = 0;
params.constantConcen=0;
params.coshQ = 1;
params.rhoConcenQ = 1;
params.concen_step = 10;

%% Concentration parameters
params.h_time_c = 0.05;        % time step size for concentration
params.t_final_c =300;        % final time for concentration

params.mu_cluster =0;%0.045;%0.045;  %0.05; %0.1 for everything except chem var_rad
params.receptor_scalar =0.01;%0.005;%0.2;%0.005;%.015;


%% Dependent parameter calculations
params.cell_radi = [params.cell_radi_nurse, params.cell_radi_oocyte, params.cell_radi_cluster];
params.x_centers = [params.x_centers_nurse, params.x_centers_oocyte, params.x_centers_cluster];
params.y_centers = [params.y_centers_nurse, params.y_centers_oocyte, params.y_centers_cluster];
params.target_volumes = [params.target_volumes_nurse, params.target_volumes_oocyte, params.target_volumes_cluster];

dof.n_space = params.x_length/params.h_space+1;  % number of points in the x and y axis
dof.n_time = params.t_final/params.h_time;       % number of time steps
dof.n_cells_nurse = length(params.cell_radi_nurse);
dof.n_cells_oocyte = length(params.cell_radi_oocyte);
dof.n_cells_cluster = length(params.cell_radi_cluster);
dof.n_cells = dof.n_cells_nurse + dof.n_cells_oocyte + dof.n_cells_cluster;

fd.x = 0:params.h_space:params.x_length;  % point vlaues along the x axis
fd.y = 0:params.h_space:params.y_length;  % point values along the y axis
[fd.x_mesh, fd.y_mesh] = meshgrid(fd.x, fd.y);  % meshgrid for the domain
fd.epithelial = tanh(((fd.y_mesh-2.5).^2 + ((fd.x_mesh-2.5)*0.65).^2)./0.5).^100;
fd.epithelial = reshape(fd.epithelial, [dof.n_space^2,1]);
fd.h_epithelial = 3*fd.epithelial.^2-2*fd.epithelial.^3;
fd.id_mat = speye(dof.n_space);
fd.laplacian_1d = (1/params.h_space^2)*spdiags(ones(dof.n_space, 1) * [1, 1, -2, 1, 1], [1-dof.n_space, -1:1, dof.n_space-1], dof.n_space, dof.n_space);
fd.laplacian = kron(fd.laplacian_1d, fd.id_mat) + kron(fd.id_mat, fd.laplacian_1d);


params.folder = fileparts(which(mfilename));    % Get the folder path the script is in
addpath(genpath(params.folder));                % Add folder plus all subfolders to the path.

%% Generate parameter string for file identification
% Create filename string based on parameters
param_string = sprintf('nr%s_x%s_y%s_v%s_o%.4f_%.1f_%.1f_%.2f_c%.4f_%.2f_%.1f_%.1f_d%s_a%.0f_bs%.1f_es%.3f_g%.2f_xl%.0f_yl%.0f_hs%.2f_ht%.2f_tf%.0f', ...
    strjoin(arrayfun(@(x) sprintf('%.3g',x), params.cell_radi_nurse, 'UniformOutput', false),'-'), ...
    strjoin(arrayfun(@(x) sprintf('%.2g',x), params.x_centers_nurse, 'UniformOutput', false),'-'), ...
    strjoin(arrayfun(@(x) sprintf('%.2g',x), params.y_centers_nurse, 'UniformOutput', false),'-'), ...
    strjoin(arrayfun(@(x) sprintf('%.2g',x), params.target_volumes_nurse, 'UniformOutput', false),'-'), ...
    params.cell_radi_oocyte, params.x_centers_oocyte, params.y_centers_oocyte, params.target_volumes_oocyte, ...
    params.cell_radi_cluster, params.x_centers_cluster, params.y_centers_cluster, params.target_volumes_cluster, ...
    strjoin(arrayfun(@(x) sprintf('%.2g',x), params.dval, 'UniformOutput', false),'-'), ...
    params.alpha, params.beta_s, params.eta_s, params.gamma, ...
    params.x_length, params.y_length, params.h_space, params.h_time, params.t_final);

% Clean up the string to be filename-safe
param_string = regexprep(param_string, '[^\w\-]', '_');
params.param_id_string = param_string;

% Clean up the string to be filename-safe
params.param_id_string = strrep(param_string, '.', 'p');
params.param_id_string = strrep(params.param_id_string, '-', 'n');

return