% EE688: Optimal Control Theory
% Fall 2022, KAIST
% Author: Hany Hamed
% Assignment 3
% Problem1: Receding Horizon Control
% Code adapted from the lecture code

clc;
clear all;
close all;

% System dynamics
A = [1.5431 1.1752; 1.1752 1.5431]; B = [0.5431 1.1752]';

xmax = 1;
X_poly = Polyhedron('lb', [-xmax;-xmax], 'ub', [xmax;xmax]);
U_poly = Polyhedron('lb', -0.5, 'ub', 0.5);

% set C constraint Ac * x <= bc
Ac = [1 1; -1 -1];
bc = [0.4;0.4];
C_poly = Polyhedron(Ac, bc);


%% (a)

% The set C is control invariant if it is subset of Pre(C)

% Precursor of C
set_poly = C_poly;
A_om = set_poly.A; 
b_om = set_poly.b;
Au = U_poly.A;
bu = U_poly.b;

Pre_CU = Polyhedron('A', [A_om * A, A_om * B ; zeros(2,2), Au], 'b', [b_om; bu]);
Pre_C = Pre_CU.projection([1 2]);


figure(1);
grid on; hold on;
plot(Pre_C, 'color', 'm');
plot(set_poly, 'color', 'r');
legend("Pre(C)", "", "", "", "C");


%% (b) N-step backward reachable set K_N(X_f=C)
close all;


figure(2)
grid on; hold on;
xlabel('x_1');
ylabel('x_2'); 
set(gca, 'fontsize', 12)
plot([-xmax xmax xmax -xmax -xmax], [-xmax -xmax xmax xmax -xmax], 'k', 'LineWidth', 2)

% set_poly = X_poly;
set_poly = C_poly;
Omega = set_poly;
A_om = set_poly.A; 
b_om = set_poly.b;
Au = U_poly.A;
bu = U_poly.b;

for i=1:100
    % Precursor of Xf=C
    Pre_XU = Polyhedron('A', [A_om * A, A_om * B ; zeros(2,2), Au], 'b', [b_om; bu]);
    Pre_X = Pre_XU.projection([1 2]);

    A_pre = Pre_X.H(:,1:2);
    b_pre = Pre_X.H(:,3);

    Omega_next = Polyhedron('A', [A_pre; A_om], 'b', [b_pre; b_om]);
    
    if Omega_next == Omega
        disp(i)
        disp('converged')
        break
    end
    
    Omega = Omega_next;
    A_om = Omega.H(:,1:2);
    b_om = Omega.H(:,3);
end

Nbar = i;

% Cinf supposdly to be the X_0 feasible region
Cinf = Omega;
% TODO: Something is wrong!!!!
plot(C_poly, 'color', 'r')
plot(Cinf, 'color', 'm')

% 
% 
% figure(2)
% grid on; hold on;
% xlabel('x_1');
% ylabel('x_2'); 
% set(gca, 'fontsize', 12)
% plot([-xmax xmax xmax -xmax -xmax], [-xmax -xmax xmax xmax -xmax], 'k', 'LineWidth', 2)
% 
% set_poly = C_poly;
% Omega = set_poly;
% A_om = set_poly.A; 
% b_om = set_poly.b;
% Au = U_poly.A;
% bu = U_poly.b;
% 
% for i=1:100
%     % Precursor of Xf=C
%     Pre_XU = Polyhedron('A', [A_om * A, A_om * B ; zeros(2,2), Au], 'b', [b_om; bu]);
%     Pre_X = Pre_XU.projection([1 2]);
% 
%     A_pre = Pre_X.H(:,1:2);
%     b_pre = Pre_X.H(:,3);
% 
%     Omega_next = Polyhedron('A', [A_pre; A_om], 'b', [b_pre; b_om]);
%     
%     if Omega_next == Omega
%         disp(i)
%         disp('converged')
%         break
%     end
%     
%     Omega = Omega_next;
%     A_om = Omega.H(:,1:2);
%     b_om = Omega.H(:,3);
% end
% 
% Nbar = i;
% Cinf = Omega;
% % plot(Pre_X, 'color', 'r')
% plot(Cinf, 'color', 'm')
%% (c)


%% Testing
A = [1.5431 1.1752;1.1752 1.5431];
B = [0.5431;1.1752];

% state constraint Ax * x <= bx
xmax = 0.4;
Ax = [1 1; -1 -1];
bx = xmax * ones(2,1);

% input constraint Au * u <= bu
umax = 0.5;
Au = [1; -1]; 
bu = umax * ones(2,1);

X = Polyhedron(Ax, bx);
U = Polyhedron(Au, bu);

% state constraint
H_X = X.H(:,1:2); K_X = X.H(:,3);

% input constraint
H_U = U.H(:,1); K_U = U.H(:,2);

% terminal state set
Xf = X;

H_Xf = Xf.H(:,1:2); K_Xf = Xf.H(:,3);
N = 1;

figure(3)
subplot(121)
hold on; axis equal;
plot(X, 'color',[0.7 0.7 0.7]);
axis([-3 3 -3 3])
xlabel('x_1');
ylabel('x_2'); 
set(gca, 'fontsize', 12)

for k=1:N
    
    % Precursor of X
    Pre_XU = Polyhedron('A', [H_Xf * A, H_Xf * B ; zeros(2,2), H_U], 'b', [K_Xf; K_U]);
    Pre_X = Pre_XU.projection([1 2]);

    % One-step BRS of X
    H_pre = Pre_X.H(:,1:2);
    K_pre = Pre_X.H(:,3);

    BRS_X = Polyhedron('A', [H_pre; H_X], 'b', [K_pre; K_X]);
    
    H_Xf = BRS_X.H(:,1:2);
    K_Xf = BRS_X.H(:,3);
    
    n = [ '1', num2str(N+1),num2str(k+1)];
    subplot(n)
    hold on; axis equal;
    
    plot(Pre_X, 'color', 'r')
%     plot(X,'color',[0.7 0.7 0.7]);

    axis([-3 3 -3 3])
    xlabel('x_1');
    ylabel('x_2'); 
    set(gca, 'fontsize', 12)

end

