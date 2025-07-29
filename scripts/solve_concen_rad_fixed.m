function concen = solve_concen_rad_fixed(x, y, u_cells, epithelial, n_space)

    t = [linspace(0, 100000, 50)];
    borders_cell =  epithelial.^2;

    for j = 1:8
        borders_cell = borders_cell + u_cells(:,j).^2;
    end
    borders_cell = ones(size(borders_cell))-borders_cell;

    brd_cell = reshape(borders_cell,n_space,n_space);

    x_dom = x(25:59);
    y_dom = y(41:61);

    sol = pdepe(0,  @(x, t, c, dcdx) pde_fun(x, t, c, dcdx), @icfun, @(xL, cL, xR, cR, t) bc_fun(xL, cL, xR, cR, t), x_dom, t);
    concen = sol(end,:);
end
function [g, f, s] = pde_fun(x, t, c, dcdx)
    % Parameters
    D = 10;
    k = 16.8;%0.00042;
    A = 1;
    g = A;
    f = D * A .* dcdx;
    s = -k *A .* c;
end


function c0 = icfun(~)
    c0 = 15;
end

function [pL, qL, pR, qR] = bc_fun(xL, cL, xR, cR, t) 
    D=10;
    % Set boundary conditions
    pL = 0;
    qL = 1;
    %A = pi * r(x_dom, radius_values, xR).^2;
    pR = -35.54;%-10./(D*A); 
    qR = 1;
end