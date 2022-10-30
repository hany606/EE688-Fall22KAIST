% Source adapted from the lecture codes
% Author: Hany Hamed
% Assignment 2: problem 3: EE688, KAIST, Fall 2022
% Linear quadraric regulator
clc;
clear all;
close all;

%% (a)

A_c = [0 1; 1 0]; B_c = [0; 1];
dt = 1;

% Analytically
%A = expm(A_c*dt); B = inv(A_c)*(A_c-eye(2))*B_c;

% Using matlab cmd
sys = ss(A_c, B_c, [1 0 ;0 1], [0;0]);
sys_d = c2d(sys, dt, 'zoh');

A = sys_d.A; B = sys_d.B;

%% (b) MPC, N=1
P = diag([1,0.1]); R = 10^3; N = 1;
x0 = [0, -0.3]';
N = 1;
% state constraint Ax * x <= bx
xmax = 1;
Ax = [eye(2); -eye(2)];
bx = xmax * ones(4,1);

% input constraint Au * u <= bu
umax = 0.5;
Au = [1; -1]; 
bu = umax * ones(2,1);

simulation_time = 100;

x = zeros(size(A,2), simulation_time+1);
u = zeros(size(B,2), simulation_time);
x(:,1) = x0;

figure(1)
grid on; hold on;
xlabel('x_1');
ylabel('x_2'); 
set(gca, 'fontsize', 12)
plot([-xmax xmax xmax -xmax -xmax], [-xmax -xmax xmax xmax -xmax], 'k', 'LineWidth', 2)


for k=1:simulation_time
    % compute the MPC control
    if k == 1
        [xstar, ustar, Jstar, exitflag ] = p3_finite_control1(P, R, A, B, Ax, bx, Au, bu, x0, [],[], 3,1);
        if exitflag ~= 1
            error('The problem is infeasible with the initial state!')
        end
    else
        [xstar, ustar, Jstar, exitflag ] = p3_finite_control1(P, R, A, B, Ax, bx, Au, bu, x(:,k), [], [], 3,1);
    end
    
    if exitflag ~= 1
        mssg = ['Infeasible at time step ', num2str(k)];
        display(mssg)
        break
    end

    x(:,k+1) = A*x(:,k) + B*ustar(1);
    
    figure(1)
    plot(x(1,1:k+1), x(2,1:k+1), '-o', 'LineWidth', 2, 'Color', 'b')
end



%% (c) Control-invariant

P = diag([1,0.1]); R = 10^3; N = 1;
x_0 = [0, -0.3]';
N = 1;

% state constraint Ax * x <= bx
Ax1 = [1 1; -1 -1];
bx1 = 0.4 * ones(2,1);


simulation_time = 20;

x = zeros(size(A,2), simulation_time+1);
x(:,1) = x0;

figure(2)
grid on; hold on;
xlabel('x_1');
ylabel('x_2'); 
set(gca, 'fontsize', 12)
plot([-xmax xmax xmax -xmax -xmax], [-xmax -xmax xmax xmax -xmax], 'k', 'LineWidth', 2)


for k=1:simulation_time
    % compute the MPC control
    if k == 1
        [xstar, ustar, Jstar, exitflag ] = p3_finite_control1(P, R, A, B, Ax, bx, Au, bu, x0, Ax1,bx1, 3,1);
        if exitflag ~= 1
            error('The problem is infeasible with the initial state!')
        end
    else
        [xstar, ustar, Jstar, exitflag ] = p3_finite_control1(P, R, A, B, Ax, bx, Au, bu, x(:,k),Ax1,bx1, 3,1);
    end
    
    if exitflag ~= 1
        mssg = ['Infeasible at time step ', num2str(k)];
        display(mssg)
        break
    end
    x(:,k+1) = A*x(:,k) + B*ustar(1);
    
    figure(2)

    plot(x(1,1:k+1), x(2,1:k+1), '-o', 'LineWidth', 2, 'Color', 'b')
end



%% (d) Feasible regions
% figure(3)
% grid on; hold on;
% xlabel('x_1');
% ylabel('x_2'); 
% set(gca, 'fontsize', 12)
% 
% feasible   = [];
% 
% feasible = [feasible; [0,3]];
% feasible = [feasible; [-1,2]];
% feasible = [feasible; [1,3]];
% 
% % feasible(:,:,2)= [-1,2];
% % 
% % 
% % feasible(:,:,3)= [1,3];
% 
% 
% disp(feasible(:,1));
% 
% 
% plot(feasible(:,1), feasible(:,2), 'ro')

N_d = 15;

simulation_time = 5;
x1_grid = linspace(-xmax,xmax,N_d);
x2_grid = linspace(-xmax,xmax,N_d);

feasible   = [];
feasible_c = [];
feasible_both = [];

figure(3)
grid on; hold on;
xlabel('x_1');
ylabel('x_2'); 
set(gca, 'fontsize', 12)
plot([-xmax xmax xmax -xmax -xmax], [-xmax -xmax xmax xmax -xmax], 'k', 'LineWidth', 2)

for i=1:N_d
    x1 = x1_grid(i);
    for j=1:N_d
        x2 = x2_grid(j);
        disp([i,j]);
        x0 = [x1, x2];
        flag1 = 0;
        x = zeros(size(A,2), simulation_time+1);
        x(:,1) = x0;
        for k=1:simulation_time
            % compute the MPC control
            if k == 1
                [xstar, ustar, Jstar, exitflag ] = p3_finite_control1(P, R, A, B, Ax, bx, Au, bu, x(:,1), [],[], 0,1);
                if exitflag ~= 1
%                     error('The problem is infeasible with the initial state!')
                    mssg = ['Infeasible at initial step ', num2str(k)];
                    display(mssg)
                    flag1 = 1;
                    break
                end
            else
                [xstar, ustar, Jstar, exitflag ] = p3_finite_control1(P, R, A, B, Ax, bx, Au, bu, x(:,k), [], [],0, 1);
            end
            
            if exitflag ~= 1
                flag1 = 1;
                mssg = ['Infeasible at time step ', num2str(k)];
                display(mssg)
                break
            end
            x(:,k+1) = A*x(:,k) + B*ustar(1);

%             figure(1)
%             plot(x(1,1:k+1), x(2,1:k+1), '-o', 'LineWidth', 2, 'Color', 'b')
        end
        

        flag2 = 0;
        x = zeros(size(A,2), simulation_time+1);
        x(:,1) = x0;
        for k=1:simulation_time
            % compute the MPC control
            if k == 1
                [xstar, ustar, Jstar, exitflag ] = p3_finite_control1(P, R, A, B, Ax, bx, Au, bu, x(:,1), Ax1,bx1, 0,1);
                if exitflag ~= 1
%                     error('The problem is infeasible with the initial state!')
                    mssg = ['Infeasible at initial step ', num2str(k)];
                    display(mssg)
                    flag2 = 1;
                    break
                end
            else
                [xstar, ustar, Jstar, exitflag ] = p3_finite_control1(P, R, A, B, Ax, bx, Au, bu, x(:,k),Ax1,bx1, 0,1);
            end
            
            if exitflag ~= 1
                flag2 = 1;
                mssg = ['Infeasible at time step ', num2str(k)];
                display(mssg)
                break
            end
            x(:,k+1) = A*x(:,k) + B*ustar(1);

%             figure(1)
%             plot(x(1,1:k+1), x(2,1:k+1), '-o', 'LineWidth', 2, 'Color', 'b')
        end
        if flag1==0 && flag2==0
            feasible_both = [feasible_both; x0];
        elseif flag1 == 0 && flag2 == 1
            feasible = [feasible; x0];
        elseif flag2 == 0 && flag1 == 1
            feasible_c = [feasible_c; x0];
        end 
    end
end
if(~isempty(feasible_both))
    plot(feasible_both(:,1), feasible_both(:,2), 'go', 'MarkerFaceColor', 'g')
end
if(~isempty(feasible_c))
    plot(feasible_c(:,1), feasible_c(:,2), 'bo', 'MarkerFaceColor', 'b')
end
if(~isempty(feasible))
    plot(feasible(:,1), feasible(:,2), 'ro', 'MarkerFaceColor', 'r')
end