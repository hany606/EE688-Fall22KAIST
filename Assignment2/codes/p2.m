% Source adapted from the lecture codes
% Author: Hany Hamed
% Assignment 2: problem 2: EE688, KAIST, Fall 2022
% Linear quadraric regulator





clc;
clear all;
close all;

% System dynamics
zeta = 0.1; omega_n = 2; delta_t = 0.1;
A = [1 delta_t; -omega_n^2*delta_t 1-2*zeta*omega_n*delta_t]; B = [0 delta_t]';
x0 = [1 1]';

simulation_time_sec = 10; % in seconds

%% LQR specific
Q = [1 0; 0 1];  R = 1;

simulation_time = simulation_time_sec/delta_t;


%% Infinite-horizon LQR (c)
P_inf = dare(A, B, Q, R);
F_inf = -inv(B'*P_inf*B+R)*B'*P_inf*A;


x = zeros(size(A,2), simulation_time+1);
u = zeros(size(B,2), simulation_time);
x(:,1) = x0;

for i=1:simulation_time
    u(:,i) = F_inf * x(:,i);
    x(:,i+1) = (A + B*F_inf) * x(:,i);
end



figure(1)
subplot(311)
grid on; hold on;
plot(x(1,:), '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('x_1');
title('Infinite horizon');
xlim([1, simulation_time+1])
set(gca, 'fontsize', 12)

subplot(312)
grid on; hold on;
plot(x(2,:), '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('x_2')
xlim([1, simulation_time+1])
set(gca, 'fontsize', 12)

subplot(313)
grid on; hold on;
plot(u, '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('u')
xlim([1, simulation_time])
set(gca, 'fontsize', 12)

%% LQR is a PD controller (b) Part1: changing Q

% Q = 5*[5 7; 2 8];  R = 1;
N_d = 100;

coeff_q = linspace(0.1,50,N_d);
lambda1 = []; lambda2 = [];
fp = []; fd=[];
area1 = [];
area2 = [];
for j=1:N_d
    q = coeff_q(j);
    Q = q * diag([1,1]); R =1;
    lambda = eig(Q);
    lambda1 = [lambda1; lambda(1)]; lambda2 = [lambda2; lambda(2)];

    P_inf = dare(A, B, Q, R);
    F_inf = -inv(B'*P_inf*B+R)*B'*P_inf*A;
    fp = [fp; F_inf(1)]; fd = [fd; F_inf(2)];

    x = zeros(size(A,2), simulation_time+1);
    u = zeros(size(B,2), simulation_time);
    x(:,1) = x0;
    
    area1_ = 0;
    area2_ = 0;
    for i=1:simulation_time
        u(:,i) = F_inf * x(:,i);
        x(:,i+1) = (A + B*F_inf) * x(:,i);
        area1_ = area1_ + abs(x(1,i+1));
        area2_ = area2_ + abs(x(2,i+1));

    end
    area1 = [area1; area1_];
    area2 = [area2; area2_];
end

%clf(2);

figure(2)
subplot(211)
grid on; hold on;
plot(coeff_q, fp, 'LineWidth', 2)
plot(coeff_q, fd, 'LineWidth', 2)
plot(coeff_q, lambda1, 'LineWidth', 2)
plot(coeff_q, lambda2, 'LineWidth', 2)
legend("F_p", "F_d", "\lambda_1", "\lambda_2")
ylabel('value');
xlabel('Q-coeffiecient');
title('Q, Fp, Fd characterstics');
subplot(212)
grid on; hold on;
plot(coeff_q, area1, 'LineWidth', 2)
plot(coeff_q, area2, 'LineWidth', 2)
legend("x_1", "x_2")
ylabel('Cumulated sum of area');
xlabel('Q-coeffiecient');
title('Cumulated sum of area under curve');

%% (b) continue plotting changed Q
% Q = 0.1 * diag([1,1]); R =1;
Q = 50 * diag([1,1]); R =1;

P_inf = dare(A, B, Q, R);
F_inf = -inv(B'*P_inf*B+R)*B'*P_inf*A;


x = zeros(size(A,2), simulation_time+1);
u = zeros(size(B,2), simulation_time);
x(:,1) = x0;

for i=1:simulation_time
    u(:,i) = F_inf * x(:,i);
    x(:,i+1) = (A + B*F_inf) * x(:,i);
end



figure(3)
subplot(311)
grid on; hold on;
plot(x(1,:), '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('x_1');
% title('Infinite horizon');
xlim([1, simulation_time+1])
set(gca, 'fontsize', 12)

subplot(312)
grid on; hold on;
plot(x(2,:), '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('x_2')
xlim([1, simulation_time+1])
set(gca, 'fontsize', 12)

subplot(313)
grid on; hold on;
plot(u, '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('u')
xlim([1, simulation_time])
set(gca, 'fontsize', 12)


%% LQR is a PD controller (b) Part2: changing R

N_d = 100;

coeff_r = linspace(0.1,50,N_d);
fp = []; fd=[];
area1 = [];
area2 = [];
for j=1:N_d
    r = coeff_r(j);
    Q = diag([1,1]); R =r;
    P_inf = dare(A, B, Q, R);
    F_inf = -inv(B'*P_inf*B+R)*B'*P_inf*A;
    fp = [fp; F_inf(1)]; fd = [fd; F_inf(2)];

    x = zeros(size(A,2), simulation_time+1);
    u = zeros(size(B,2), simulation_time);
    x(:,1) = x0;
    
    area1_ = 0;
    area2_ = 0;
    for i=1:simulation_time
        u(:,i) = F_inf * x(:,i);
        x(:,i+1) = (A + B*F_inf) * x(:,i);
        area1_ = area1_ + abs(x(1,i+1));
        area2_ = area2_ + abs(x(2,i+1));

    end
    area1 = [area1; area1_];
    area2 = [area2; area2_];
end

%clf(2);

figure(4)
subplot(211)
grid on; hold on;
plot(coeff_r, fp, 'LineWidth', 2)
plot(coeff_r, fd, 'LineWidth', 2)
legend("F_p", "F_d")
ylabel('value');
xlabel('R');
title('Fp, Fd characterstics');
subplot(212)
grid on; hold on;
plot(coeff_r, area1, 'LineWidth', 2)
plot(coeff_r, area2, 'LineWidth', 2)
legend("x_1", "x_2")
ylabel('Cumulated sum of area');
xlabel('Q-coeffiecient');
title('Cumulated sum of area under curve');
%% (b) continue plotting changed R
% Q = diag([1,1]); R =0.1;
Q = diag([1,1]); R =50;

P_inf = dare(A, B, Q, R);
F_inf = -inv(B'*P_inf*B+R)*B'*P_inf*A;


x = zeros(size(A,2), simulation_time+1);
u = zeros(size(B,2), simulation_time);
x(:,1) = x0;

for i=1:simulation_time
    u(:,i) = F_inf * x(:,i);
    x(:,i+1) = (A + B*F_inf) * x(:,i);
end



figure(5)
subplot(311)
grid on; hold on;
plot(x(1,:), '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('x_1');
% title('Infinite horizon');
xlim([1, simulation_time+1])
set(gca, 'fontsize', 12)

subplot(312)
grid on; hold on;
plot(x(2,:), '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('x_2')
xlim([1, simulation_time+1])
set(gca, 'fontsize', 12)

subplot(313)
grid on; hold on;
plot(u, '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('u')
xlim([1, simulation_time])
set(gca, 'fontsize', 12)
%% (b) continue plotting changed R and Q
% Q = 0.1*diag([1,1]); R =50;
Q = 50*diag([1,1]); R =0.1;

P_inf = dare(A, B, Q, R);
F_inf = -inv(B'*P_inf*B+R)*B'*P_inf*A;


x = zeros(size(A,2), simulation_time+1);
u = zeros(size(B,2), simulation_time);
x(:,1) = x0;

for i=1:simulation_time
    u(:,i) = F_inf * x(:,i);
    x(:,i+1) = (A + B*F_inf) * x(:,i);
end



figure(6)
subplot(311)
grid on; hold on;
plot(x(1,:), '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('x_1');
% title('Infinite horizon');
xlim([1, simulation_time+1])
set(gca, 'fontsize', 12)

subplot(312)
grid on; hold on;
plot(x(2,:), '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('x_2')
xlim([1, simulation_time+1])
set(gca, 'fontsize', 12)

subplot(313)
grid on; hold on;
plot(u, '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('u')
xlim([1, simulation_time])
set(gca, 'fontsize', 12)

%% Receding-horizon (c)
Q = [1 0; 0 1];  R = 1;

x = zeros(size(A,2), simulation_time+1);
u = zeros(size(B,2), simulation_time);
x(:,1) = x0;
% clf(3);
N = 30;
for i=1:simulation_time
    % compute the finite horizon
    P_N = Q;
    
    P = zeros(2,2,N+1);
    P(:,:,N+1) = P_N;
    F = zeros(size(B,2), size(B,1), N);
    for j=1:N
        k = N - j + 1;
        F(:,:,k) = -inv(B'*P(:,:,k+1)*B+R)*B'*P(:,:,k+1)*A;
        P(:,:,k) = A'* P(:,:,k+1) * A + Q - A' * P(:,:,k+1) * B * inv(B'*P(:,:,k+1)*B + R) * B' *P(:,:,k+1) * A;
    end

    u(:,i) = F(:,:,1) * x(:,i);
    x(:,i+1) = (A + B*F(:,:,1)) * x(:,i);
end

%clf(3)
figure(7)
subplot(311)
grid on; hold on;
plot(x(1,:), '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('x_1');
title('Receding-Horizon');
xlim([1, simulation_time+1])
set(gca, 'fontsize', 12)

subplot(312)
grid on; hold on;
plot(x(2,:), '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('x_2')
xlim([1, simulation_time+1])
set(gca, 'fontsize', 12)

subplot(313)
grid on; hold on;
plot(u, '-o', 'LineWidth', 2, 'Color', 'b','MarkerSize',3)
ylabel('u')
xlim([1, simulation_time])
set(gca, 'fontsize', 12)
