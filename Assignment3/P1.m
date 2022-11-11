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
Ax = X_poly.A;
bx = X_poly.b;
U_poly = Polyhedron('lb', -0.5, 'ub', 0.5);

% set C constraint Ac * x <= bc
Ac = [1 1; -1 -1];
bc = [0.4;0.4];
C_poly = Polyhedron([Ac; Ax], [bc; bx]);


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

% One-step BRS of X
A_pre = Pre_C.H(:,1:2);
b_pre = Pre_C.H(:,3);

intersection = Polyhedron([A_pre; A_om], [b_pre; b_om]); % intersection of the precursor and the original


figure(1);
grid on; hold on;
plot(Pre_C, 'color', 'r');
plot(C_poly, 'color', [0.7 0.7 0.7]);
plot(intersection, 'color', 'm'); % intersection between C and Pre(C)
legend("Pre(C)", "C", "$C \cap Pre(C)$",'Interpreter','latex');
xlabel('x_1');
ylabel('x_2'); 
set(gca, 'fontsize', 12)
fprintf("Is C in Pre(C) (Pre(C).contains(C))? -> %d\n", Pre_C.contains(C_poly)); % 1


% P1 = Polyhedron('A',[0,1;0,-1],'b',[1,1]);  % Unbounded
% P2 = Polyhedron('A',[0,1;0,-1; -1,0],'b',[1,1,1]); % Bounded from the
% left
% P1.plot('color', 'r')
% P2.plot('color', 'g')
% P1.contains([-2;0]) % 1
% P2.contains([-2;0]) % 0

%% (b) N-step backward reachable set K_N(X_f=C); C is control invariant
% close all;


figure(2)
grid on; hold on;
xlabel('x_1');
ylabel('x_2'); 
set(gca, 'fontsize', 12)
plot([-xmax xmax xmax -xmax -xmax], [-xmax -xmax xmax xmax -xmax], 'k', 'LineWidth', 2)

% set_poly = X_poly;
set_poly = C_poly;
Omega = set_poly;
A_om = Omega.A; 
b_om = Omega.b;
Au = U_poly.A;
bu = U_poly.b;

N = 1;
for i=1:N
    % Precursor of Xf=C
    Pre_XU = Polyhedron('A', [A_om * A, A_om * B ; zeros(2,2), Au], 'b', [b_om; bu]);
    Pre_X = Pre_XU.projection([1 2]);

    % One-step BRS of X
    A_pre = Pre_X.H(:,1:2);
    b_pre = Pre_X.H(:,3);

    BRS_X = intersect(X_poly, Pre_X);%Polyhedron('A', [A_pre; Ax], 'b', [b_pre; bx]);
    
    A_om = BRS_X.H(:,1:2);
    b_om = BRS_X.H(:,3);
end

plot(Pre_X, 'color', 'r');
plot(BRS_X, 'color', 'm');
plot(C_poly, 'color', [0.7 0.7 0.7]);

% X_0 is a superset of C

fprintf("Is X0 in C (C.contains(X0))? -> %d\n", C_poly.contains(BRS_X)); % 0
fprintf("Is C in X0 (X0.contains(C))? -> %d\n", BRS_X.contains(C_poly)); % 1
fprintf("Is X0 in Pre(C) (Pre(C).contains(X0))? -> %d\n", Pre_X.contains(BRS_X)); % 1

legend("$\mathcal{X}$ constraints", "Pre(C)",  "$X_0=K_1(C)=BRS_1(C)=Pre(C)\cap \mathcal{X}$", "C",'Interpreter','latex');
%% (c)


% Compute the maximal invariant set C_inf
% Omega = X_poly;
% A_om = Ax; 
% b_om = bx;
% 
% for i=1:100
%     % Precursor of X
%     Pre_XU = Polyhedron('A', [A_om * A, A_om * B ; zeros(2,2), Au], 'b', [b_om; bu]);
%     Pre_X = Pre_XU.projection([1 2]);
% 
%     A_pre = Pre_X.H(:,1:2);
%     b_pre = Pre_X.H(:,3);
% 
%     Omega_next = Polyhedron('A', [A_pre; A_om], 'b', [b_pre; b_om]);
%     
%     if Omega_next == Omega
%         disp('converged at')
%         disp(i);
%         break
%     end
%     
%     Omega = Omega_next;
%     A_om = Omega.H(:,1:2);
%     b_om = Omega.H(:,3);
% end
% Cinf = Omega;

Q = diag([1,0.1]); R=1000;
P_infty= dare(A, B, Q, R);
P = P_infty;

F_infty = -inv(B'*P_infty*B+R)*B'*P_infty*A;

A_auto = A + B * F_infty;

figure(3)
grid on; hold on;
xlabel('x_1');
ylabel('x_2'); 
set(gca, 'fontsize', 12)
plot([-xmax xmax xmax -xmax -xmax], [-xmax -xmax xmax xmax -xmax], 'k', 'LineWidth', 2)

Omega = Polyhedron([Ax; Au * F_infty], [bx; bu]);

A_om = [Ax; Au * F_infty]; 
b_om = [bx; bu];

for i=1:100
    % Precursor of X
    Pre_X = Polyhedron('A', [A_om * A_auto; Au * F_infty], 'b', [b_om; bu]);

    A_pre = Pre_X.H(:,1:2);
    b_pre = Pre_X.H(:,3);

    Omega_next = Polyhedron('A', [A_pre; A_om], 'b', [b_pre; b_om]);
    
    if Omega_next == Omega
        disp('converged at')
        disp(i);
        break
    end
    
    Omega = Omega_next;
    A_om = Omega.H(:,1:2);
    b_om = Omega.H(:,3);
end

Nbar = i;
O_infty = Omega;
% plot(Cinf, 'color', 'm')
plot(O_infty, 'color', 'y')
% legend("X constraints", "$C_\infty$", "$O_\infty$",'Interpreter','latex');
legend("X constraints", "$O_\infty$",'Interpreter','latex');


simulation_time = 20;
x = zeros(size(A,2), simulation_time+1);
x0 = [0, -0.3]'; % From assignment 2
x(:,1) = x0;
Axf = O_infty.A;
bxf = O_infty.b;

figure(4)
grid on; hold on;
xlabel('x_1');
ylabel('x_2'); 
set(gca, 'fontsize', 12)
plot([-xmax xmax xmax -xmax -xmax], [-xmax -xmax xmax xmax -xmax], 'k', 'LineWidth', 2)


for k=1:simulation_time
    % compute the MPC control
    if k == 1
        [xstar, ustar, Jstar, exitflag ] = Extra.finite_control1(P, R, A, B, Ax, bx, Au, bu, x0, Axf,bxf, 3,1);
        if exitflag ~= 1
            error('The problem is infeasible with the initial state!')
        end
    else
        [xstar, ustar, Jstar, exitflag ] = Extra.finite_control1(P, R, A, B, Ax, bx, Au, bu, x(:,k),Axf,bxf, 3,1);
    end
    
    if exitflag ~= 1
        mssg = ['Infeasible at time step ', num2str(k)];
        display(mssg)
        break
    end
    x(:,k+1) = A*x(:,k) + B*ustar(1);
    
    figure(4)

    plot(x(1,1:k+1), x(2,1:k+1), '-o', 'LineWidth', 2, 'Color', 'b')
end