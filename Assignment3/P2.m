% EE688: Optimal Control Theory
% Fall 2022, KAIST
% Author: Hany Hamed
% Tube-based MPC

clc;
clear all;
close all;

%% (a)
% Determine tightened state and input constraint sets X¯ and U¯ s
% Nominal system xk+1 = ¯xk + ¯uk satisfies the tightened

clc;
clear all;
close all;

A = 1; B=1;
K = -0.5;

% x-constraints: <= 2
X_poly = Polyhedron('ub',2);
% disp(X_poly.contains(-1000));

% u-constraints: -1 <= u <= 1
U_poly = Polyhedron('lb', -1, 'ub', 1);


% figure(1);
% hold("on");
% plot(X_poly,'color','r');
% plot(U_poly,'color','b');


% W as set
% w-constraints: -0.1 <= w <= 0.1
W_poly = Polyhedron('lb', -0.1, 'ub', 0.1);

S_hat_infty = W_poly;
for i=1:100
    % It is possible to directly use + and -
    S_hat_infty = S_hat_infty.plus(((A+B*K)^i)*W_poly);
end

figure(2);
hold("on");
plot(S_hat_infty, 'color', 'm');
plot(W_poly, 'color', 'k');

X_bar_poly = X_poly.minus(S_hat_infty);
% TODO: Is it correct that U_bar did not change?
U_bar_poly = U_poly.minus(K*S_hat_infty);

figure(3);
hold("on");
plot(X_bar_poly, 'color', 'r');
plot(U_bar_poly, 'color', 'b');

% plot(S_hat_infty, 'color', 'r');
% plot(U_poly- K*S_hat_infty, 'color', 'b');

%% (b)
% Determine tightened state and input constraint sets X¯ and U¯ s
% Nominal system xk+1 = ¯xk + ¯uk satisfies the tightened

clc;
clear all;
close all;

A = [1, 1; 0, 1]; B=[0;1];
K = [-0.4, -1.2];

% x-constraints:
X_poly = Polyhedron('lb', [-10;-10], 'ub', [10;10]);
% disp(X_poly.contains(-1000));

% u-constraints:
U_poly = Polyhedron('lb', -1, 'ub', 1);


figure(1);
hold("on");
plot(X_poly,'color','r');
plot(U_poly,'color','b');


% W as set
% w-constraints: -0.1 <= w <= 0.1
W_poly = Polyhedron('lb', [-0.1; -0.1], 'ub', [0.1; 0.1]);

S_hat_infty = W_poly;
for i=1:100
    S_hat_infty = S_hat_infty.plus(((A+B*K)^i)*W_poly);
end

figure(2);
hold("on");
plot(S_hat_infty, 'color', 'm');
plot(W_poly, 'color', 'k');

X_bar_poly = X_poly.minus(S_hat_infty);
U_bar_poly = U_poly.minus(K*S_hat_infty);

figure(3);
hold("on");
plot(X_bar_poly, 'color', 'r');
plot(U_bar_poly, 'color', 'b');



% Question: What is meant with that : mpt3 toolbox to construct a tube S^_N; N=3 ????????????????????? TODO


S_hat_N = W_poly;
N = 3;
for i=1:N
    S_hat_N = S_hat_N.plus(((A+B*K)^i)*W_poly);
end


X_bar_poly = X_poly.minus(S_hat_N);
U_bar_poly = U_poly.minus(K*S_hat_N);

figure(4);
hold("on");
plot(X_bar_poly, 'color', 'r');
plot(U_bar_poly, 'color', 'b');

%% (c)
% Question: TODO: what is the initial condition
% Question: TODO: Is there any terminal constrain I should care about?

% close all;

% Nominal system
% MPC for the system

P = eye(2); Q = eye(2); R=1;
Ax = X_bar_poly.A;
bx = X_bar_poly.b;


Au = U_bar_poly.A;
bu = U_bar_poly.b;

Af = Ax;
bf = bx;
x0 = [7;-2];



N_predicition_horizon = 3;
simulation_time = 20;
x_oinf = zeros(size(x0,1), simulation_time+1);
x_open = zeros(size(x0,1), N+1, simulation_time);
ustar_seq = zeros(1, simulation_time);
x_oinf(:,1) = x0;
figure(6);
hold("on");
for k=1:simulation_time
    w = zeros(2,1, N);
    % Solve the finite horizon optimal control problem
    if k == 1
        [xstar, ustar, Jstar, exitflag ] = Extra.finite_control( P, Q, R, A, B, Ax, bx, Au, bu, Af, bf, N_predicition_horizon, x0, w)
        if exitflag ~= 1
            error('The problem is infeasible with the initial state!')
        end
    else
        [xstar, ustar, Jstar, exitflag ] = Extra.finite_control( P, Q, R, A, B, Ax, bx, Au, bu, Af, bf, N_predicition_horizon, x_oinf(:,k), w)
    end
    
    if exitflag ~= 1
        mssg = ['Infeasible at time step ', num2str(k)];
        display(mssg)
        break
    end
    
    x_open(:,:,k) = xstar;
    ustar_seq(k) = ustar(1);
    x_oinf(:,k+1) = A*x_oinf(:,k) + B*ustar(1);
    

    plot(x_oinf(1,1:k+1), x_oinf(2,1:k+1), '-o', 'LineWidth', 2, 'Color', 'g')
    tube = Extra.reduce_poly(S_hat_N + x_oinf(:,k));
    fill(tube.V(:, 1), tube.V(:, 2),  'r', 'FaceAlpha', 0.3);
%     plot(x_open(1,:,k), x_open(2,:,k), 'or--', 'LineWidth', 1.5)

end


%% (d)

% Question: TODO: should I sample the random uncertainity again inside MPC
% and when evolve the system??

simulation_time = 20;
x_oinf1 = zeros(size(x0,1), simulation_time+1);
x_open = zeros(size(x0,1), N+1, simulation_time);
x_oinf1(:,1) = x0;
figure(7);
hold("on");
for k=1:simulation_time
    w = unifrnd(-0.1, 0.1, 2,1,3);
    % Solve the finite horizon optimal control problem
    if k == 1
        [xstar, ustar, Jstar, exitflag ] = Extra.finite_control( P, Q, R, A, B, Ax, bx, Au, bu, Af, bf, N_predicition_horizon, x0, w)
        if exitflag ~= 1
            error('The problem is infeasible with the initial state!')
        end
    else
        [xstar, ustar, Jstar, exitflag ] = Extra.finite_control( P, Q, R, A, B, Ax, bx, Au, bu, Af, bf, N_predicition_horizon, x_oinf1(:,k), w)
    end
    
    if exitflag ~= 1
        mssg = ['Infeasible at time step ', num2str(k)];
        display(mssg)
        break
    end
    
    x_open(:,:,k) = xstar;
    ustar_new = ustar_seq(k) + K*(x_oinf1(:,k) - x_oinf(:,k));
    x_oinf1(:,k+1) = A*x_oinf1(:,k) + B*ustar_new(1) + unifrnd(-0.1, 0.1, 2,1,1);
    

    plot(x_oinf1(1,1:k+1), x_oinf1(2,1:k+1), '-o', 'LineWidth', 2, 'Color', 'g')
    tube = Extra.reduce_poly(S_hat_N + x_oinf(:,k));
    fill(tube.V(:, 1), tube.V(:, 2),  'r', 'FaceAlpha', 0.3);
%     plot(x_open(1,:,k), x_open(2,:,k), 'or--', 'LineWidth', 1.5)

end

%% (d)
figure(8);
hold("on");
plot(x_oinf1(1:k+1)-x_oinf(1:k+1))


