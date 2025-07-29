% GO_SOLVER: Main solver script for phase field model simulation
%
% USAGE:
%   Set parameters in config.m
%   Visualize results using go_show_results.m
%   data files are saved in the results/ directory
%   Visualization files are saved in the visuals/ directory
%
% DEPENDENCIES:
%   - config.m: Configuration parameters
%   - solve_pfm.m: Solves for cell state without concentration
%   - solve_pfm_concen.m: Solves for cell state with concentration
%                           (requires steady state from solve_pfm.m)
%
% By Naghmeh Akhavan March 2025

clear variables

[fd, params, dof] = config();

if params.run_steadyQ
    tstart = tic;
    [sol_cells_saved, residuals] = solve_pfm(fd, params, dof);
    cd results
    save([params.param_id_string, '.mat']);
    cd ..
    fprintf('Time to solve the problem to steady state: %f\n', toc(tstart));
end

cd results
if exist([params.param_id_string, '.mat'], 'file')
    load([params.param_id_string, '.mat']);
else
    cd ..
    error('No results file found. \n Run steady state solver first via changing params.run_steadyQ to 1.');
end
cd ..
[fd, params, dof] = config();
if params.run_concenQ
    tstart = tic;
    cell_state = sol_cells_saved(:,:,:,end-1);
    [sol_cells_saved_c, residuals_c, concen, diff_coef, receptor, chemo,tension_x,tension_y] = solve_pfm_concen(fd, params, dof, cell_state);
    cd results
    save([params.param_id_string, 'concen.mat']);
    cd ..
    fprintf('Time to solve the problem to concentration: %f\n', toc(tstart));
end

return
