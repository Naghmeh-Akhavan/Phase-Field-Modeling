%% USER PARAMETERS
use_concentration =0;
type_visual =2;
type_f =3;
showTensionVecQ = 1;
resol =1; %resolution
framerate = 20;

[fd, params, dof] = config();


% Visualization options:
% type_vis = 1: 3D video
% type_vis = 2: Contour bounds video
% type_vis = 3: Contourf video
% type_vis = 4: Residuals plot   ( this option isIndependent of type_f)
% type_vis = 5: Cluster dynamics (average x positions and speeds)
%               ( this option is independent of type_f )

% What to visualize:
% type_f = 1: Cells
% type_f = 2: Concentration
% type_f = 3: Interfacial tension (TIM)
% type_f = 4: Chemical Force (F_chem)

% Concentration options:
% use_concentration = 0: No concentration
% use_concentration = 1: With concentration

cd results

% Time disretization parameters based on concentration
if use_concentration
    % results with concentration
    load([params.param_id_string, 'concen.mat']);
    total_step = dof.n_saved_c;                 % With concentration
    t_final = params.t_final_c;
    step_time = t_final/total_step;
else
    load([params.param_id_string, '.mat']);      % No concentration
    total_step = dof.n_saved;                    % No concentration
    t_final = params.t_final;
    step_time = t_final/total_step;
    show_fun = sol_cells_saved;
end

cd ..
cd visuals

% Choose the appropriate function to visualize based on type_f
switch type_f
    case 1
        if use_concentration
            show_fun = sol_cells_saved_c;
        else
            show_fun = sol_cells_saved;
        end
        v_name = 'cells';
    case 2
        show_fun = concen;
        v_name = 'concentration';
        if ~use_concentration
           error('controdicting parameters: use_concentration == 1 and type_f == 2 (concentration) is chosen');
        end
    case 3
        show_fun = receptor;
        v_name = 'interfacial_tension';
    case 4
        show_fun = chemo;
        v_name = 'chemo';
end

% Set up limits for visualization
clim_min = min(show_fun,[],'all');
clim_max = max(show_fun,[],'all');
if clim_min == clim_max
    clim_max = clim_min + 1;
    fprintf('Warning: Chosen function is constant, setting default limits.\n Consider changing the type_f.\n');
end

switch type_visual
    case 1
        % Create a VideoWriter object
        v = VideoWriter([v_name '_3D.avi']);
        v.FrameRate = framerate;
        open(v);

        figure('units','pixels','position',[0 0 560*resol 420*resol]);
        hold on
        switch type_f
            case 1
                f_epi = surf(fd.x, fd.y, reshape(fd.epithelial,dof.n_space,dof.n_space), 'FaceAlpha', 0.7);
                f_n1 = surf(fd.x, fd.y, squeeze(show_fun(:,:,1,1)), 'FaceAlpha', 0.7);
                f_n2 = surf(fd.x, fd.y, squeeze(show_fun(:,:,2,1)), 'FaceAlpha', 0.7);
                f_n3 = surf(fd.x, fd.y, squeeze(show_fun(:,:,3,1)), 'FaceAlpha', 0.7);
                f_n4 = surf(fd.x, fd.y, squeeze(show_fun(:,:,4,1)), 'FaceAlpha', 0.7);
                f_n5 = surf(fd.x, fd.y, squeeze(show_fun(:,:,5,1)), 'FaceAlpha', 0.7);
                f_n6 = surf(fd.x, fd.y, squeeze(show_fun(:,:,6,1)), 'FaceAlpha', 0.7);
                f_oo = surf(fd.x, fd.y, squeeze(show_fun(:,:,7,1)), 'FaceAlpha', 0.7);
                f_cl = surf(fd.x, fd.y, squeeze(show_fun(:,:,8,1)), 'FaceAlpha', 0.7);
                title(sprintf('Time: %.2f', 0));
                shading("interp")
                view(240,60);
                axis tight;
                for i = 2:total_step
                    set(f_n1, 'ZData', squeeze(show_fun(:,:,1,i)));
                    set(f_n2, 'ZData', squeeze(show_fun(:,:,2,i)));
                    set(f_n3, 'ZData', squeeze(show_fun(:,:,3,i)));
                    set(f_n4, 'ZData', squeeze(show_fun(:,:,4,i)));
                    set(f_n5, 'ZData', squeeze(show_fun(:,:,5,i)));
                    set(f_n6, 'ZData', squeeze(show_fun(:,:,6,i)));
                    set(f_oo, 'ZData', squeeze(show_fun(:,:,7,i)));
                    set(f_cl, 'ZData', squeeze(show_fun(:,:,8,i)));
                    % Update title with current time
                    title(sprintf('Time: %.2f', i*step_time));
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end

            case 2
                f_c = surf(fd.x, fd.y, squeeze(show_fun(:,:,1)));
                title(sprintf('Concentration at Time: %.2f', 0));
                view(240,60);
                axis tight;
                for i = 2:total_step
                    set(f_c, 'ZData', squeeze(show_fun(:,:,i)));
                    % Update title with current time
                    title(sprintf('Concentration at Time: %.2f', i*step_time));
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end

            case 3
                f_r = surf(fd.x, fd.y, squeeze(show_fun(:,:,1)));
                title(sprintf('Interfacial Tension at Time: %.2f', 0));
                view(240,80);
                axis tight;
                for i = 2:total_step
                    set(f_r, 'ZData', squeeze(show_fun(:,:,i)));
                    % Update title with current time
                    title(sprintf('Interfacial Tension at Time: %.2f', i*step_time));
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end
              case 4
                f_chemo = surf(fd.x, fd.y, squeeze(show_fun(:,:,1)));
                title(sprintf('Interfacial Tension at Time: %.2f', 0));
                view(240,80);
                axis tight;
                for i = 2:total_step
                    set(f_chemo, 'ZData', squeeze(show_fun(:,:,i)));
                    % Update title with current time
                    title(sprintf('Chemical Force at Time: %.2f', i*step_time));
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end
        end
        % Close the video file
        close(v);

    case 2
        % Create a VideoWriter object
        v = VideoWriter([v_name '_cont_bnd.avi']);
        v.FrameRate = framerate;
        open(v);

        figure('units','pixels','position',[0 0 560*resol 420*resol]);
        hold on
        switch type_f
            case 1
                kk =1;
                f_epi = contour(fd.x, fd.y, reshape(fd.epithelial,dof.n_space,dof.n_space),[0.5,0.5],'black',"LineWidth",3);
                [f_n1, h_n1] = contour(fd.x, fd.y, squeeze(show_fun(:,:,1,kk)),[0.5,0.5], 'red',"LineWidth",1.5);
                [f_n2, h_n2] = contour(fd.x, fd.y, squeeze(show_fun(:,:,2,kk)),[0.5,0.5], 'red',"LineWidth",1.5);
                [f_n3, h_n3] = contour(fd.x, fd.y, squeeze(show_fun(:,:,3,kk)),[0.5,0.5], 'red',"LineWidth",1.5);
                [f_n4, h_n4] = contour(fd.x, fd.y, squeeze(show_fun(:,:,4,kk)),[0.5,0.5], 'red',"LineWidth",1.5);
                [f_n5, h_n5] = contour(fd.x, fd.y, squeeze(show_fun(:,:,5,kk)),[0.5,0.5], 'red',"LineWidth",1.5);
                [f_n6, h_n6] = contour(fd.x, fd.y, squeeze(show_fun(:,:,6,kk)),[0.5,0.5], 'red',"LineWidth",1.5);
                [f_oo, h_oo] = contour(fd.x, fd.y, squeeze(show_fun(:,:,7,kk)),[0.5,0.5],'blue',"LineWidth",1.5);
                [f_cl, h_cl] = contour(fd.x, fd.y, squeeze(show_fun(:,:,8,kk)),[0.5,0.5],'green',"LineWidth",1.5);

                title(sprintf('Time: %.2f', kk*step_time));
                axis tight;
                if showTensionVecQ
                    q1 = quiver(fd.x, fd.y, tension_x(:,:,1), tension_y(:,:,1));
                end
                for i = 2:total_step
                    set(h_n1, 'ZData', squeeze(show_fun(:,:,1,i)));
                    set(h_n2, 'ZData', squeeze(show_fun(:,:,2,i)));
                    set(h_n3, 'ZData', squeeze(show_fun(:,:,3,i)));
                    set(h_n4, 'ZData', squeeze(show_fun(:,:,4,i)));
                    set(h_n5, 'ZData', squeeze(show_fun(:,:,5,i)));
                    set(h_n6, 'ZData', squeeze(show_fun(:,:,6,i)));
                    set(h_oo, 'ZData', squeeze(show_fun(:,:,7,i)));
                    set(h_cl, 'ZData', squeeze(show_fun(:,:,8,i)));
                    % Update title with current time
                    title(sprintf('Time: %.2f', i*step_time));
                    if showTensionVecQ
                        set(q1, 'UData', tension_x(:,:,i), 'VData', tension_y(:,:,i));
                    end
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end

            case 2
                [f_c, h_c] = contour(fd.x, fd.y, squeeze(show_fun(:,:,1)));
                clim([clim_min clim_max]);
                title(sprintf('Concentration at Time: %.2f', 0));
                axis tight;
                for i = 2:total_step
                    set(h_c, 'ZData', squeeze(show_fun(:,:,i)));
                    % Update title with current time
                    title(sprintf('Concentration at Time: %.2f', i*step_time));
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end

            case 3
                [f_r, h_r] = contour(fd.x, fd.y, squeeze(show_fun(:,:,1)));
                clim([clim_min clim_max]);
                title(sprintf('Interfacial Tension at Time: %.2f', 0));
                axis tight;
                for i = 2:total_step
                    set(h_r, 'ZData', squeeze(show_fun(:,:,i)));
                    % Update title with current time
                    title(sprintf('Interfacial Tension at Time: %.2f', i*step_time));
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end
            case 4
                [f_chemo, h_chemo] = contour(fd.x, fd.y, squeeze(show_fun(:,:,1)));
                clim([clim_min clim_max]);
                title(sprintf('Chemical Force at Time: %.2f', 0));
                axis tight;
                for i = 2:total_step
                    set(h_chemo, 'ZData', squeeze(show_fun(:,:,i)));
                    % Update title with current time
                    title(sprintf('Chemical Force at Time: %.2f', i*step_time));
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end
        end
        % Close the video file
        close(v);

    case 3
        % Create a VideoWriter object
        v = VideoWriter([v_name '_contourf.avi']);
        v.FrameRate = framerate;
        open(v);

        figure('units','pixels','position',[0 0 560*resol 420*resol]);
        hold on
        switch type_f
            case 1
                f_epi = contour(fd.x, fd.y, reshape(fd.epithelial,dof.n_space,dof.n_space),[0.5,0.5],'black',"LineWidth",3);
                [f_n1, h_n1] = contourf(fd.x, fd.y, squeeze(show_fun(:,:,1,1)),  0.5:0.1:1,LineStyle="none");
                [f_n2, h_n2] = contourf(fd.x, fd.y, squeeze(show_fun(:,:,2,1)),  0.5:0.1:1,LineStyle="none");
                [f_n3, h_n3] = contourf(fd.x, fd.y, squeeze(show_fun(:,:,3,1)),  0.5:0.1:1,LineStyle="none");
                [f_n4, h_n4] = contourf(fd.x, fd.y, squeeze(show_fun(:,:,4,1)),  0.5:0.1:1,LineStyle="none");
                [f_n5, h_n5] = contourf(fd.x, fd.y, squeeze(show_fun(:,:,5,1)),  0.5:0.1:1,LineStyle="none");
                [f_n6, h_n6] = contourf(fd.x, fd.y, squeeze(show_fun(:,:,6,1)),  0.5:0.1:1,LineStyle="none");
                [f_oo, h_oo] = contourf(fd.x, fd.y, squeeze(show_fun(:,:,7,1)),  0.5:0.1:1,LineStyle="none");
                [f_cl, h_cl] = contourf(fd.x, fd.y, squeeze(show_fun(:,:,8,1)),  0.5:0.1:1,LineStyle="none");
                title(sprintf('Time: %.2f', 0));
                colormap(flipud(gray));
                shading("interp")
                axis tight;
                for i = 2:total_step
                    set(h_n1, 'ZData', squeeze(show_fun(:,:,1,i)));
                    set(h_n2, 'ZData', squeeze(show_fun(:,:,2,i)));
                    set(h_n3, 'ZData', squeeze(show_fun(:,:,3,i)));
                    set(h_n4, 'ZData', squeeze(show_fun(:,:,4,i)));
                    set(h_n5, 'ZData', squeeze(show_fun(:,:,5,i)));
                    set(h_n6, 'ZData', squeeze(show_fun(:,:,6,i)));
                    set(h_oo, 'ZData', squeeze(show_fun(:,:,7,i)));
                    set(h_cl, 'ZData', squeeze(show_fun(:,:,8,i)));
                    % Update title with current time
                    title(sprintf('Time: %.2f', i*step_time));
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end

            case 2
                [f_c, h_c] = contourf(fd.x, fd.y, squeeze(show_fun(:,:,80)));
                clim([clim_min clim_max]);
                title(sprintf('Concentration at Time: %.2f', 80*step_time));
                axis tight;
                for i = 2:total_step
                    set(h_c, 'ZData', squeeze(show_fun(:,:,i)));
                    % Update title with current time
                    title(sprintf('Concentration at Time: %.2f', i*step_time));
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end

            case 3
                [f_r, h_r] = contourf(fd.x, fd.y, squeeze(show_fun(:,:,1)));
                clim([clim_min clim_max]);
                title(sprintf('Interfacial Tension at Time: %.2f', 0));
                axis tight;
                for i = 2:total_step
                    set(h_r, 'ZData', squeeze(show_fun(:,:,i)));
                    % Update title with current time
                    title(sprintf('Interfacial Tension at Time: %.2f', i*step_time));
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end
            case 4
                [f_chemo, h_chemo] = contourf(fd.x, fd.y, squeeze(show_fun(:,:,1)));
                clim([clim_min clim_max]);
                title(sprintf('Chemical Force at Time: %.2f', 0));
                axis tight;
                for i = 2:total_step
                    set(h_chemo, 'ZData', squeeze(show_fun(:,:,i)));
                    % Update title with current time
                    title(sprintf('Chemical Force at Time: %.2f', i*step_time));
                    frame = getframe(gcf);
                    writeVideo(v, frame);
                end
        end
        % Close the video file
        close(v);
    case 4
        % Plot residuals
        figure('units','pixels','position',[0 0 560*resol 420*resol]);
        plot(1:total_step, residuals_c, 'LineWidth', 2);
        xlabel('Time Step');
        ylabel('Residuals');
        title('Residuals over Time Steps');
        grid on;
    case 5
        % Calculate weighted average x positions and speeds for cluster cells over time
        u_cells_cluster = sol_cells_saved_c(:, :, 8, :);
        avg_x_positions = zeros(1, total_step-1);
        speeds = zeros(1, total_step-2);
        for t = 1:(total_step-1)
            total_weight = 0;
            weighted_sum = 0;
            for i = 1:dof.n_space
                weight = sum(u_cells_cluster(:, i, t));
                weighted_sum = weighted_sum + weight * fd.x(i);
                total_weight = total_weight + weight;
            end
            avg_x_positions(t) = weighted_sum / total_weight;
            if t > 1
                speeds(t-1) = (avg_x_positions(t) - avg_x_positions(t-1)) / step_time;
            end
        end

        % Create time vector
        time = (1:(total_step-1))*step_time;

        % Plot average X position
        figure(1);
        plot(time, avg_x_positions, 'LineWidth', 3);
        ylabel('\bfAverage X Position', 'FontWeight', 'bold', 'FontSize', 14);
        xlabel('\bfTime', 'FontWeight', 'bold', 'FontSize', 14);
        title('\bfWeighted Average X Position of Cluster Cell Over Time', 'FontWeight', 'bold', 'FontSize', 14);
        grid off;
        saveas(gcf, 'avg_x_position_cluster.jpg', 'jpg');

        % Plot speed
        figure(2);
        plot(time(1:end-1), speeds, 'LineWidth', 3);
        ylabel('\bfSpeed', 'FontWeight', 'bold', 'FontSize', 14);
        xlabel('\bfTime', 'FontWeight', 'bold', 'FontSize', 14);
        title('\bfSpeed of Cluster Cell Over Time', 'FontWeight', 'bold', 'FontSize', 14);
        grid off;
        saveas(gcf, 'speed_cluster.jpg', 'jpg');

        % Plot average X position and speed on the same figure with separate y-axes
        figure(3);
        yyaxis left
        plot(time, avg_x_positions, 'b-', 'LineWidth', 3);
        ylabel('\bfAverage X Position', 'Color', 'b', 'FontWeight', 'bold', 'FontSize', 14);
        xlabel('\bfTime', 'FontWeight', 'bold', 'FontSize', 14);
        title('\bfWeighted Average X Position and Speed of Cluster Cell Over Time', 'FontWeight', 'bold', 'FontSize', 14);
        set(gca, 'YColor', 'b');

        yyaxis right
        plot(time(1:end-1), speeds, 'r-', 'LineWidth', 3);
        ylabel('\bfSpeed', 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 14);
        set(gca, 'YColor', 'r');

        grid on;
        legend('Average X Position', 'Speed', 'Location', 'best');
        saveas(gcf, 'position_speed_cluster.jpg', 'jpg');
  

end

cd ..