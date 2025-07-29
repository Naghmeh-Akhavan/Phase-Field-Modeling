function concen = solve_concen_rad(x, y, u_cells, epithelial, n_space)

    t = [linspace(0, 100000, 50)];
    borders_cell =  epithelial.^2;
            
    for j = 1:8
        borders_cell = borders_cell + u_cells(:,j).^2;
    end
    borders_cell = ones(size(borders_cell))-borders_cell;

    brd_cell = reshape(borders_cell,n_space,n_space);

    x_dom = x(25:59);
    y_dom = y(41:61);
    radius_values = r_vals(x_dom, y_dom, brd_cell);

    sol = pdepe(0,  @(x, t, c, dcdx) pde_fun(x, t, c, dcdx, x_dom, radius_values), @icfun, @(xL, cL, xR, cR, t) bc_fun(xL, cL, xR, cR, t, x_dom, radius_values), x_dom, t);
    concen = sol(end,:);
end
function [g, f, s] = pde_fun(x, t, c, dcdx, x_dom, radius_values)
    % Parameters
    D = 14;      
    k = 16.8;%0.00042; 
    A = pi *  r(x_dom, radius_values, x).^2; % Cross-sectional area
    g = A;
    f = D * A .* dcdx;
    s = -k *A .* c;
end


function c0 = icfun(~)
    c0 = 15;
end

function [pL, qL, pR, qR] = bc_fun(xL, cL, xR, cR, t, x_dom, radius_values) 
    D=1;
    % Set boundary conditions
    pL = 0; 
    qL = 1;
    %A = pi * r(x_dom, radius_values, xR).^2;
    pR = -1;%-10./(D*A); 
    qR = 1; 
end

function rad_val = r(x_dom, raduis_values, x_val)
    rad_val = interp1(x_dom, raduis_values, x_val,'spline','extrap');
end

function radius_values = r_vals(x, y_dom, brd_cell)
    %r = rb * ones(size(x));
    n_x = length(x);
    x_1 = zeros(n_x,1);
    x_2 = zeros(n_x,1);
    target_val = 0.6;
    for j = 1:n_x
        temp = brd_cell(41:61,j+25);
        [maxValue, index] = max(temp);
        for i = 1:(index-1)
            if x_1(j) == 0 && temp(index - i) < target_val
                x_1(j) = y_dom(index - i);
            end
        end
        for i = 1:(21-index)
            if x_2(j) == 0 && temp(index + i) < target_val
                x_2(j) = y_dom(index + i);
            end
        end
    end
    radius_values = (x_2 - x_1)./2;
end
