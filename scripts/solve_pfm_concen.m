function [sol_cells_saved, residuals, concen_saved,diff_coef_saved, receptor_saved,chemo_saved, tension_x_saved,tension_y_saved] = solve_pfm_concen(fd, params, dof, cell_state)
% Solve the PFM model with concentration
% Inputs:
%   fd:         struct containing the mesh and finite diffrence information
%   params:     struct containing the parameters of the model
%   dof:        struct containing the degrees of freedom
%   cell_state: saved solutions of the PFM model without concentration
%
% Outputs:
%   sol_cells_saved:    saved solutions of the PFM model
%   residuals:          residuals of the PFM model
%
% By Naghmeh Akhavan March 2025

dof.n_saved = dof.n_saved_c;
params.h_time = params.h_time_c;
params.t_final = params.t_final_c;
dof.n_time = params.t_final/params.h_time; % Number of time steps fixed for concentration run

% Solution to the cells (:,:,1:n) nurse cells, (:,:,end-1) oocyte, (:,:,end) cluster
sol_cells = reshape(cell_state,dof.n_space^2,8);

sol_cells_saved = zeros(dof.n_space,dof.n_space, dof.n_cells, dof.n_saved);
concen_saved = zeros(dof.n_space,dof.n_space, dof.n_saved);
diff_coef_saved = zeros(dof.n_space,dof.n_space, dof.n_saved);
receptor_saved = zeros(dof.n_space,dof.n_space, dof.n_saved);
chemo_saved = zeros(dof.n_space,dof.n_space, dof.n_saved);
tension_x_saved =  zeros(dof.n_space,dof.n_space, dof.n_saved);
tension_y_saved = zeros(dof.n_space,dof.n_space, dof.n_saved);
residuals = [];

% Central difference with Neumann boundary conditions
id_mat = speye(dof.n_space);

central_diff_1d = -(1/(2*params.h_space))*spdiags(ones(dof.n_space, 1) * [ 1, -1], [ -1,1], dof.n_space, dof.n_space);
central_diff_1d(1,1) = -1/(2*params.h_space);
central_diff_1d(dof.n_space,dof.n_space) = 1/(2*params.h_space);
central_diff_x = kron(central_diff_1d, id_mat);
central_diff_y = kron(id_mat, central_diff_1d);

% diffusion coefficient
mu_cluster =params.mu_cluster;
saved_iter = 0;
percent = 0;
saved_t = tic;

if params.fixed_ConcenQ
    concen = solve_concen_rad_fixed(fd.x, fd.y, sol_cells, fd.epithelial, dof.n_space);
else
    concen = solve_concen_rad(fd.x, fd.y, sol_cells, fd.epithelial, dof.n_space);
end

concen = concen.^3./((concen.^2+0.01).*(concen+0.05));
concen= kron([concen(1)*ones(24,1);concen';concen(end)*ones(dof.n_space-59,1)],ones(dof.n_space,1));
concen = reshape(concen, dof.n_space, dof.n_space);
[concen_x, concen_y] = gradient(concen,params.h_space);
concen_x = reshape(concen_x, [dof.n_space^2, 1]);
concen_y = reshape(concen_y, [dof.n_space^2, 1]);
concen = reshape(concen, dof.n_space^2, 1);


for i = 1:dof.n_time
    sol_cells = round(sol_cells, 16);
    sol_cells_old = sol_cells;
    % Compute h(u)= 3u^2-2u^3
    h_sol_cells = 3*sol_cells.^2-2*sol_cells.^3;

    % phi(u) = \sum h(u)
    phi_sum = sum(h_sol_cells(:,1:end-2),2);
    vol_cells = params.h_space^2*sum(h_sol_cells);

    % Nurse cells Evaluation
    for m = 1:dof.n_cells_nurse
        sol_cells(:,m) = sol_cells(:,m) + params.h_time*(...
            params.dval(1)*fd.laplacian*sol_cells(:,m)...
            + sol_cells(:,m).*(1-sol_cells(:,m)).*(...
                sol_cells(:,m)-0.5+(params.target_volumes(m)-vol_cells(m))...
                -params.beta(1, 1)*(phi_sum - h_sol_cells(:,m))...   %NC-NC
                -params.beta(1, 2)*h_sol_cells(:,end-1)...                              %NC-Oocyte
                -params.beta(1, 3)*h_sol_cells(:,end)...                             %NC-Clus
            - params.beta_s*fd.h_epithelial...
                +fd.laplacian*(...
            + params.eta_s*fd.h_epithelial...
            + params.eta(1, 1)*(phi_sum - h_sol_cells(:,m))...
            + params.eta(1, 2)*h_sol_cells(:,end-1)...
            + params.eta(1, 3)*h_sol_cells(:,end)...
            + params.gamma*h_sol_cells(:,m))));
    end

    % Oocyte Evaluation
    sol_cells(:,end-1) =  sol_cells(:,end-1) + params.h_time*(...
        params.dval(2)*fd.laplacian* sol_cells(:,end-1)...
        +  sol_cells(:,end-1).*(1- sol_cells(:,end-1)).*(...
            sol_cells(:,end-1)-0.5+(params.target_volumes(end-1)-vol_cells(end-1))...
            -params.beta(2, 1)*phi_sum...                               %oocy-NC
            -params.beta(2, 3)*h_sol_cells(:,end)...                             %oocy-clus
        - params.beta_s*fd.h_epithelial...
            +fd.laplacian*(...
        + params.eta_s*fd.h_epithelial...
        + params.eta(2, 1)*phi_sum...
        + params.eta(2, 3)*h_sol_cells(:,end)...
        + params.gamma*h_sol_cells(:,end-1))));

    if params.run_tensionQ
        [cluster_x, cluster_y] = gradient(reshape(sol_cells(:,end),dof.n_space,dof.n_space),params.h_space);
        cluster_x = reshape(cluster_x, [dof.n_space^2, 1]);
        cluster_y = reshape(cluster_y, [dof.n_space^2, 1]);
        if params.coshQ
            rho_concen = concen;
        else
            rho_concen =concen.^3./((concen.^2+1).*(concen+1));
        end
        tension_x = (rho_concen.*sol_cells(:,end).*sum(sol_cells(:,1:6),2).*(sign(-concen_y.*cluster_x+concen_x.*cluster_y).*cluster_y)) ;
        tension_y = -(rho_concen.*sol_cells(:,end).*sum(sol_cells(:,1:6),2).*(sign(-concen_y.*cluster_x+concen_x.*cluster_y).*cluster_x));
        receptorF = -params.receptor_scalar*(...
            +central_diff_x*(tension_x) ...
            +central_diff_y*(tension_y) ...
        );
    else
        receptorF = 0;
    end
    % Cluster Evaluation
    chem_x = sol_cells(:,end).*concen_x;
    chem_y = sol_cells(:,end).*concen_y;
    sol_cells(:,end) =  sol_cells(:,end) + params.h_time*(...
        +receptorF...
        -mu_cluster*(central_diff_x*chem_x + central_diff_y*(chem_y))...
        +params.dval(3)*fd.laplacian* sol_cells(:,end)...
        +  sol_cells(:,end).*(1- sol_cells(:,end)).*(...
            sol_cells(:,end)-0.5+params.alpha*(params.target_volumes(end)-vol_cells(end))...
            -params.beta(3, 1)*phi_sum...                               %clus-NC
            -params.beta(3, 2)*h_sol_cells(:,end-1)...                             %clus-oocy
        - params.beta_s*fd.h_epithelial... 
            +fd.laplacian*(...
        + params.eta_s*fd.h_epithelial...
        + params.eta(3, 1)*phi_sum...
        + params.eta(3, 2)*h_sol_cells(:,end-1)...
        + params.gamma*h_sol_cells(:,end))));

    res = norm(sol_cells(:,end) - sol_cells_old(:,end))/norm(sol_cells_old(:,end));
    if mod(i, dof.n_time/dof.n_saved) == 0
        saved_iter = saved_iter + 1;
        sol_cells_saved(:,:,:,i/(dof.n_time/dof.n_saved)) = reshape(sol_cells, dof.n_space, dof.n_space, dof.n_cells);
        concen_saved(:,:,i/(dof.n_time/dof.n_saved)) = reshape(concen, dof.n_space, dof.n_space);
        if params.coshQ == 0
           diff_coef_saved(:,:,i/(dof.n_time/dof.n_saved)) = reshape(diff_coef, dof.n_space, dof.n_space);
        end
        if params.run_tensionQ
            receptor_saved(:,:,i/(dof.n_time/dof.n_saved)) = reshape(receptorF, dof.n_space, dof.n_space);
            tension_x_saved(:,:,i/(dof.n_time/dof.n_saved))= reshape(tension_x, dof.n_space, dof.n_space);
            tension_y_saved(:,:,i/(dof.n_time/dof.n_saved))= reshape(tension_y, dof.n_space, dof.n_space);
        else
            tension_x_saved(:,:,i/(dof.n_time/dof.n_saved))= reshape(chem_x, dof.n_space, dof.n_space);
            tension_y_saved(:,:,i/(dof.n_time/dof.n_saved))= reshape(chem_y, dof.n_space, dof.n_space);
        end

        chemo_saved(:,:,i/(dof.n_time/dof.n_saved)) = reshape(-mu_cluster*(central_diff_x*(sol_cells(:,end).*concen_x) + central_diff_y*(sol_cells(:,end).*concen_y)), dof.n_space, dof.n_space);
        fprintf('Time step %d of %d is saved. Relative residual for cluster %d Tip of cluster %d \n', i, dof.n_time, res, max(sol_cells_old(:,end),[],'all'));
        if saved_iter == 50
            saved_iter = 0;
            percent = percent +10;
            toc(saved_t);
            fprintf('Percent %d is completed. Estimated time remaining %d seconds \n', percent, (toc(saved_t)) * ((100-percent)/10));
            saved_t = tic;
        end
        residuals = [residuals, res];
        if res < params.tol
            break
        end
    end
end

return