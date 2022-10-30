% The code is adapted from the lecture codes to include only a single
% step MPC
function [ xstar, ustar, Jstar, exitflag ] = p3_finite_control1( P, R, A, B, Ax, bx, Au, bu, x0, Ax1,bx1, verbose, debug)

    yalmip('clear')
    N = 1;
    % setup the problem
    u = sdpvar(1,N);
    x = sdpvar(2,N+1);

    x(:,1) = x0;
    
    objective = 0;
    k = 1;
    constraints =  [Au* u(k)<= bu, Ax * (A*x(:,k)+B*u(k))<= bx];
    if(~isempty(Ax1))
        constraints = [constraints, Ax1 * (A*x(:,k)+B*u(k))<= bx1];
    end
    objective = objective + (A*x(:,k)+B*u(k))'*P*(A*x(:,k)+B*u(k)) + u(k)'*R*u(k);

    % solve the problem
    ops = sdpsettings;
    ops.solver = 'quadprog';
    ops.verbose = verbose;
    ops.debug = debug;
    diagnostics = optimize(constraints, objective, ops)

    % results
    if diagnostics.problem == 0
        exitflag= 1;
        ustar = value(u);
        xstar = value(x);
        Jstar = value(objective);
    else
        exitflag = 0;
        ustar = [];
        xstar = [];
        Jstar = [];
    end

end