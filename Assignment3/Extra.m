% EE688: Optimal Control Theory
% Fall 2022, KAIST
% Author: Hany Hamed
% Assignment 3
% Extra codes for utilities

classdef Extra    
    methods(Static)
        
        function reduced = reduce_poly(poly)
            vert = poly.V;
            x_vert = round(vert(:, 1), 5);
            y_vert = round(vert(:, 2), 5);
            idx = convhull(x_vert, y_vert);
            reduced = Polyhedron([x_vert(idx), y_vert(idx)]);
        end

        function [ xstar, ustar, Jstar, exitflag ] = finite_control( P, Q, R, A, B, Ax, bx, Au, bu, Af, bf, N, x0, w)
            yalmip('clear')
            % setup the problem
            u = sdpvar(1,N);
            x = sdpvar(2,N+1);
        
            x(:,1) = x0;
            if ~all(Ax*x(:,1)<= bx + 1e-9)
                exitflag = 0;
                ustar = [];
                xstar = [];
                Jstar = [];
                return
            end
            objective = 0;
            constraints = [];
            for k = 1:N
                if k == 1
                    constraints = [constraints, x(:,k+1) == A*x(:,k) + B*u(k) + w(:,k), Au* u(k)<= bu];
                else
                    constraints = [constraints, x(:,k+1) == A*x(:,k) + B*u(k) + w(:,k), Au* u(k)<= bu, Ax * x(:,k)<= bx];
                end
                objective = objective + x(:,k)'*Q*x(:,k) + u(k)'*R*u(k);
            end
            constraints = [constraints, Ax * x(:,N+1) <= bx];
            if ~isempty(Af)
                constraints = [constraints, Af * x(:,N+1) <= bf];
            end
            objective = objective + x(:,N+1)'*P*x(:,N+1);
        
            
            % solve the problem
            ops = sdpsettings;
            ops.solver = 'quadprog';
            ops.verbose = 3;
            ops.debug = 1;
        
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

        % The code is adapted from the lecture codes to include only a single
        % step MPC (From Assignment 2)
        function [ xstar, ustar, Jstar, exitflag ] = finite_control1( P, R, A, B, Ax, bx, Au, bu, x0, Ax1,bx1, verbose, debug)
        
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


    end
end

