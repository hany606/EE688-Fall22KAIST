% EE688: Optimal Control Theory
% Fall 2022, KAIST
% Author: Hany Hamed
% Assignment 3
% Problem2: Tube-based MPC

clc;
clear all;
close all;

%% (a)
% Determine tightened state and input constraint sets X¯ and U¯ s
% Nominal system xk+1 = ¯xk + ¯uk satisfies the tightened

A = 1; B=1;
K = -0.5;

% x-constraints: <= 2
X_poly = Polyhedron('ub',2);
% disp(X_poly.contains(-1000));

% u-constraints: -1 <= u <= 1
U_poly = Polyhedron('lb', -1, 'ub', 1);


figure('Name', "(a)");
grid on; hold on;
plot(X_poly,'color','r');
plot(U_poly,'color','b');
legend("", "X constraints", "U constraints");
xlabel('x_1');
ylabel('x_2'); 
title("Constraints")


% W as set
% w-constraints: -0.1 <= w <= 0.1
W_poly = Polyhedron('lb', -0.1, 'ub', 0.1);

S_hat_infty = W_poly;
S_hat_next = S_hat_infty;
% TODO: Make infinity by the convergence of the S_hat
for i=1:100
    S_hat_infty = S_hat_infty.plus(((A+B*K)^i)*W_poly);
    if(S_hat_next == S_hat_infty)
        disp("Converge at");
        disp(i);
        break
    end
    S_hat_next = S_hat_infty;
end

figure('Name', "(a)");
grid on; hold on;
plot(S_hat_infty, 'color', 'm');
plot(W_poly, 'color', 'k');
legend("$\hat{S}_\infty$", "W constraints",'Interpreter','latex');
xlabel('x_1');
ylabel('x_2'); 
title("Debug")

X_bar_poly = X_poly.minus(S_hat_infty);
% TODO: Is it correct that U_bar did not change?
U_bar_poly = U_poly.minus(K*S_hat_infty);

figure('Name', "(a)");
grid on; hold on;
plot(X_bar_poly, 'color', 'r');
plot(U_bar_poly, 'color', 'b');
legend("","$\bar{X}$", "$\bar{U}$",'Interpreter','latex');
xlabel('x_1');
ylabel('x_2'); 
title("Tightened constraints")

%% (b)
% Determine tightened state and input constraint sets X¯ and U¯ s
% Nominal system xk+1 = ¯xk + ¯uk satisfies the tightened

% clc;
% clear all;
% close all;

A = [1, 1; 0, 1]; B=[0;1];
K = [-0.4, -1.2];

% x-constraints:
X_poly = Polyhedron('lb', [-10;-10], 'ub', [10;10]);
% disp(X_poly.contains(-1000));

% u-constraints:
U_poly = Polyhedron('lb', -1, 'ub', 1);


figure('Name', "(b)");
grid on; hold on;
plot(X_poly,'color','r');
plot(U_poly,'color','b');
legend("X constraints", "U constraints");
xlabel('x_1');
ylabel('x_2'); 
title("Constraints")

% W as set
% w-constraints: -0.1 <= w <= 0.1
W_poly = Polyhedron('lb', [-0.1; -0.1], 'ub', [0.1; 0.1]);

S_hat_infty = W_poly;
S_hat_next = S_hat_infty;
for i=1:100
    S_hat_infty = S_hat_infty.plus(((A+B*K)^i)*W_poly);
    if(S_hat_next == S_hat_infty)
        disp("Converge at");
        disp(i);
        break
    end
    S_hat_next = S_hat_infty;
end

figure('Name', "(b)");
grid on; hold on;
plot(S_hat_infty, 'color', 'm');
plot(W_poly, 'color', 'k');
legend("$\hat{S}_\infty$", "W constraints",'Interpreter','latex');
xlabel('x_1');
ylabel('x_2'); 
title("Debug")


X_bar_poly = X_poly.minus(S_hat_infty);
U_bar_poly = U_poly.minus(K*S_hat_infty);

figure('Name', "(b)");
grid on; hold on;
plot(X_bar_poly, 'color', 'r');
plot(U_bar_poly, 'color', 'b');
legend("$\bar{X}$", "$\bar{U}$",'Interpreter','latex');
xlabel('x_1');
ylabel('x_2'); 
title("Tightened constraints using $\hat{S}_\infty$",'Interpreter','latex')

% Question: What is meant with that : mpt3 toolbox to construct a tube S^_N; N=3 ????????????????????? TODO


S_hat_N = W_poly;
N = 3;
for i=1:N
    S_hat_N = S_hat_N.plus(((A+B*K)^i)*W_poly);
end


figure('Name', "(b)");
grid on; hold on;
plot(S_hat_infty, 'color', 'm');
plot(W_poly, 'color', 'k');
legend("$\hat{S}_3$", "W constraints",'Interpreter','latex');
xlabel('x_1');
ylabel('x_2'); 
title("Debug")

X_bar_poly = X_poly.minus(S_hat_N);
U_bar_poly = U_poly.minus(K*S_hat_N);

figure('Name', "(b)");
grid on; hold on;
plot(X_bar_poly, 'color', 'r');
plot(U_bar_poly, 'color', 'b');
legend("$\bar{X}$", "$\bar{U}$",'Interpreter','latex');
xlabel('x_1');
ylabel('x_2'); 
title("Tightened constraints using $\hat{S}_3$",'Interpreter','latex')


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
ustar_seq = zeros(1, simulation_time);
x_oinf(:,1) = x0;
figure(9);
grid on; hold on;
xlabel('x_1');
ylabel('x_2'); 
title("Nominal system (MPC)",'Interpreter','latex')
for k=1:simulation_time
    w = zeros(2,1, N_predicition_horizon);
    % Solve the finite horizon optimal control problem
    if k == 1
        [xstar, ustar, Jstar, exitflag ] = Extra.finite_control( P, Q, R, A, B, Ax, bx, Au, bu, Af, bf, N_predicition_horizon, x0, w);
        if exitflag ~= 1
            error('The problem is infeasible with the initial state!')
        end
    else
        [xstar, ustar, Jstar, exitflag ] = Extra.finite_control( P, Q, R, A, B, Ax, bx, Au, bu, Af, bf, N_predicition_horizon, x_oinf(:,k), w);
    end
    
    if exitflag ~= 1
        mssg = ['Infeasible at time step ', num2str(k)];
        display(mssg)
        break
    end
    
    ustar_seq(k) = ustar(1);
    x_oinf(:,k+1) = A*x_oinf(:,k) + B*ustar(1);
    
    figure(9)
    plot(x_oinf(1,1:k+1), x_oinf(2,1:k+1), '-o', 'LineWidth', 2, 'Color', 'g')
    tube_plt = Extra.reduce_poly(S_hat_N + x_oinf(:,k));
    fill(tube_plt.V(:, 1), tube_plt.V(:, 2),  'r', 'FaceAlpha', 0.3);

end
legend("Nominal trajectory", "Tube");%,'Interpreter','latex');


%% (d) using the nominal trajectory from (c)
simulation_time = 20;
x_oinf1 = zeros(size(x0,1), simulation_time+1);
errors = zeros(size(x0, 1), simulation_time+1);
tube_size = zeros(2, simulation_time);
x_oinf1(:,1) = x0;
figure(10);
grid on; hold on;
xlabel('x_1');
ylabel('x_2'); 
title("Tube-based MPC",'Interpreter','latex')
for k=1:simulation_time
    % Tube-based MPC using the nominal trajectory used from (c)
    % TODO: Should we recalclate ustar_seq(k) based on the new x or from the
    % nominal trajcetory (No)
    ustar_new = ustar_seq(k) + K*(x_oinf1(:,k) - x_oinf(:,k));
    x_oinf1(:,k+1) = A*x_oinf1(:,k) + B*ustar_new(1) + unifrnd(-0.1, 0.1, 2,1,1);
    errors(:, k+1) = x_oinf1(:,k+1) - x_oinf(:,k+1);
    figure(10)

    plot(x_oinf(1,1:k+1), x_oinf(2,1:k+1), '--o', 'LineWidth', 2, 'Color', 'g')
    plot(x_oinf1(1,1:k+1), x_oinf1(2,1:k+1), '-o', 'LineWidth', 2, 'Color', 'g')
    tube_vol = S_hat_N.volume;
    tube_size(1,k) = -tube_vol/2;
    tube_size(2,k) = tube_vol/2;

    tube_plt = Extra.reduce_poly(S_hat_N + x_oinf(:,k));
    fill(tube_plt.V(:, 1), tube_plt.V(:, 2),  'r', 'FaceAlpha', 0.3);


%     tube_plt = Extra.reduce_poly(S_hat_infty + x_oinf(:,k));
%     fill(tube_plt.V(:, 1), tube_plt.V(:, 2),  'm', 'FaceAlpha', 0.3);
end
legend("Nominal trajectory", "Trajectory Tracking using Linear Feedback control", "Tube");%,'Interpreter','latex');

% N_predicition_horizon = 19;
% Question: TODO: When I make the system tracking using MPC, We need to
% reformulate the cost function to work on minimizing the error dynamics
% with the system constraints.
% should I sample the random uncertainity again inside MPC
% and when evolve the system??

% for k=1:simulation_time
%     w = unifrnd(-0.1, 0.1, 2,1,N_predicition_horizon);
%     % Solve the finite horizon optimal control problem
%     if k == 1
%         [xstar, ustar, Jstar, exitflag ] = Extra.finite_control( P, Q, R, A, B, Ax, bx, Au, bu, Af, bf, N_predicition_horizon, x0, w);
%         if exitflag ~= 1
%             error('The problem is infeasible with the initial state!')
%         end
%     else
%         [xstar, ustar, Jstar, exitflag ] = Extra.finite_control( P, Q, R, A, B, Ax, bx, Au, bu, Af, bf, N_predicition_horizon, x_oinf1(:,k), w);
%     end
%     
%     if exitflag ~= 1
%         mssg = ['Infeasible at time step ', num2str(k)];
%         display(mssg)
%         break
%     end
%     
%     x_open(:,:,k) = xstar;
%     ustar_new = ustar_seq(k) + K*(x_oinf1(:,k) - x_oinf(:,k));
%     x_oinf1(:,k+1) = A*x_oinf1(:,k) + B*ustar_new(1) + unifrnd(-0.1, 0.1, 2,1,1);
%     
%     plot(x_oinf(1,1:k+1), x_oinf(2,1:k+1), '--o', 'LineWidth', 2, 'Color', 'g')
%     plot(x_oinf1(1,1:k+1), x_oinf1(2,1:k+1), '-o', 'LineWidth', 2, 'Color', 'g')
%     tube = Extra.reduce_poly(S_hat_N + x_oinf(:,k));
%     fill(tube.V(:, 1), tube.V(:, 2),  'r', 'FaceAlpha', 0.3);
% %     plot(x_open(1,:,k), x_open(2,:,k), 'or--', 'LineWidth', 1.5)
% 
% end


%% (e)
% close all;
figure('Name', "(e)");
grid on; hold on;
% The plot of the polytope is with grey in order to show better the vectors
plot(S_hat_N, 'color', [0.7,0.7,0.7]);
scatter(errors(1,:), errors(2,:), 'color', 'b');
plotv(errors,'-o');
xlabel('x_1');
ylabel('x_2'); 
title("$\hat{S}_3$ plot with errors",'Interpreter','latex')


% figure(9);
% grid on; hold on;
% % The plot of the polytope is with grey in order to show better the vectors
% scatter(errors(1,:), errors(2,:), 'color', 'b');
% xlabel('simulation times');
% ylabel('x_2'); 
% title("Errors",'Interpreter','latex')
