%
% MPC + LSCR for ARX models
%

addpath('cop/');

%%%%% System setup. a1 and b1 are unknown. %%%%%
a1 = 1.4;
a2 = 1;
b1 = 0.6;
b2 = 1.6;
% Time
T = 100;
% Initial state
y = zeros(1,T);
y(1) = 8;
y(2) = 7;
% Time instants of kind k * T_recalc when a new set is calculated using
% LCSR_ARX
T_recalc = 10;

% Epsilon-set dimenstions. That is if the set of pissible parameters is
% rectangle with size a, b than precision is a / N, b / N.
N = 50;
% Omega is a set of the model uncertainties.
% a1 in [Omega_AB{1}(1,1), Omega(1,2)]
% b1 in [Omega(2,1), Omega(2,2)]
% w in [Omega(3,1), Omega(3,2)]
% But initially it is equal to the whole Omega.
Omega_AB = {[0 2; 0 2], N, ones(N, N)};
Omega_W = [0 5];
% C_u are constraints for the control u -- interval
C_u = [-10 10];
% C_x are constraints for the state. Here it is supposed
C_y = [-20 20];
% Name of solver to use
cop_solver = 'solve_rmpc_cop2';

fprintf('Solver used: %s\n', cop_solver);

% Noise is unknown but bounded and located within [Omega(3,1), Omega(3,2)]
mu = 2;
sigma = 3;
w = normrnd(mu,sigma,[1 T]);
w(w > Omega_W(2)) = mu;
w(w < Omega_W(1)) = mu;

%%%%% MPC setup %%%%%
N_mpc = 8; % horizon
M_mpc = 10; % number of scenarios

%%%%% LSCR_ARX Setup %%%%%
M = 100;
r = 2;

% Now calculate theta bounds
ab_bounds = Omega_AB{1};
theta_bounds = [ab_bounds(2,1) ab_bounds(2,2);
                b2 - ab_bounds(1,2)*ab_bounds(2,2) ...
                b2 - ab_bounds(1,1)*ab_bounds(2,1)];

% For plotting
[theta_x, theta_y] = get_epsilon_set(theta_bounds, N);
[a_x, b_y] = get_epsilon_set(ab_bounds, N);

%%%%% Below we do the main part - MPC + LSCR_ARX %%%%%
C = [a1 a2 b1 b2];

% Calculate re-parameterization
theta0 = b1;
theta1 = (b2 - a1 * b1);

figure('Name', 'State', 'Position', [10 500 600 290], 'Color', 'w');
figure(1);
figure('Name', 'State approximation', 'Position', [10 10 600 290], 'Color', 'w');
figure(2);
figure('Name', 'Result in terms of theta0 and theta1', 'Position', [650 500 600 290], 'Color', 'w');
figure(3);
figure('Name', 'Final result for the current step', 'Position', [650 10 600 290], 'Color', 'w');
figure(4);
figure('Name', 'State', 'Position', [650 10 600 290], 'Color', 'w');
figure(5);

% FIXME Virtual control
D = sign(randn(1,T));
D(D == 0) = 1;
D = 15 * D;

v_mpc = ones(1,1);

% Do three steps to derive initial value of x0
y(3) = -a1*y(2) - a2*y(1) + b1*6 + b2*6 + w(3);
x0 = ones(1,2);
x0(1) = y(2);
x0(2) = y(3) + a1*y(2) - b1*6;

x = zeros(2, 2);
x(:, 2) = x0';

for t=2:T
    %%%%%%%%%%%%%%%%%%% MPC step start %%%%%%%%%%%%%%%%%%%

    [v_mpc, x, res] = RMPC(str2func(cop_solver), C, C_u, C_y, Omega_AB, ...
                        Omega_W, w, v_mpc, x, N_mpc, M_mpc, t);
    if res < 0
        fprintf('Error: RMPC failed\n');
        return;
    end
    
    figure(1);
    hold on;
    plot(1:1:t, v_mpc, 'Color', 'blue'); hold on;
    plot(1:1:t, v_mpc + D(1,1:t), 'Color', 'red');
    ylabel('v_t (blue), \Delta_t (red)');
    xlabel('Time (T)');
    hold off;
    
    figure(5);
    hold on;
    plot(1:1:t, y(1,1:t), 'Color', 'blue');
    ylabel('State y_t');
    xlabel('Time (T)');
    hold off;

    %%%%%%%%%%%%%%%%%%% MPC step end %%%%%%%%%%%%%%%%%%%

    u = v_mpc;

    % Measure state FIXME : t > 2
    if t > 2
        y(t+1) = b1*u(t) + b2*u(t-1) - a1*y(t) - a2*y(t-1) + w(t+1);
    end

    x(1,t+1) = y(t+1);

    if t >= 10 && mod(t, T_recalc) == 0
        %%%%%%%%%%%%%%%%%%% LSCR step start %%%%%%%%%%%%%%%%%%%
        result0 = LSCR_ARX((t-2)/r,N,M,y,u,w,D,0,r,1,theta0,theta1,2,theta_bounds);
        result1 = LSCR_ARX((t-2)/r,N,M,y,u,w,D,1,r,1,theta0,theta1,2,theta_bounds);

        result_common = result0.*result1;

        % result3 = LSCR_ARX((t-2)/r,N,M,y,u,w,D,0,r,0,theta0,theta1,1,Omega_AB{1});

        result_ab = theta_to_ab(b2, ab_bounds, theta_bounds, N, result_common);

        Omega_AB{3} = result_ab;

        figure(3);
        hold on;
        pcolor(theta_x, theta_y, result_common);
        plot(theta0, theta1, 'r*');
        xlabel('theta0 (= b0)');
        ylabel('theta1 (= b1 - a0 * b0)');
        hold off;
        
        figure(4);
        hold on;
        pcolor(a_x, b_y, result_ab);
        plot(a1, b1, 'r*');
        xlabel('a1');
        ylabel('b1');
        hold off;
        %%%%%%%%%%%%%%%%%%% LSCR step end %%%%%%%%%%%%%%%%%%%
    end

    drawnow
end

function [x, y] = get_epsilon_set(bounds, N)
    a = bounds(1,1);
    b = bounds(1,2);
    c = bounds(2,1);
    d = bounds(2,2);
    x = a:((b-a)/N):(b-(b-a)/N);
    y = c:((d-c)/N):(d-(d-c)/N);
end

function res_ab = theta_to_ab(b2, bounds_ab,bounds_theta,N,theta_set)
% Calculate unknown parameters a1 and b1 from theta0 and theta1
    res_ab = zeros(N,N);
    a = bounds_ab(1,1);
    b = bounds_ab(1,2);
    c = bounds_ab(2,1);
    d = bounds_ab(2,2);
    a1 = bounds_theta(1,1);
    b1 = bounds_theta(1,2);
    c1 = bounds_theta(2,1);
    d1 = bounds_theta(2,2);
    for i=1:N
        for j=1:N
            a1_val = a + (b-a)*i/N;
            b1_val = c + (d-c)*j/N;

            t0 = b1_val;
            t1 = b2 - a1_val * b1_val;

            theta_i = round(N*(t0-a1)/(b1-a1));
            theta_j = round(N*(t1-c1)/(d1-c1));

            if (theta_i > 0) && (theta_i < N) ...
                    && (theta_j > 0) && (theta_j < N) ...
                    && theta_set(theta_i,theta_j) == 1
                res_ab(i,j) = 1;
            end
        end
    end
end