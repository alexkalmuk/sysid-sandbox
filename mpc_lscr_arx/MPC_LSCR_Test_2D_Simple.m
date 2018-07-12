%
% MPC + LSCR for ARX models
%

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
C_y = [-8 8];

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
theta = 0:(2/N):(2-2/N);

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

    [v_mpc, x] = RMPC(C, C_u, C_y, Omega_AB, Omega_W, w, v_mpc, x, N_mpc, M_mpc, t);
    
    figure(1);
    hold on;
    
    plot(1:1:t, v_mpc); hold on;
    plot(1:1:t, v_mpc + D(1,1:t));
    ylabel('v_t (blue), \Delta_t (red)');
    xlabel('Time (T)');
    pause(0.01);
    hold off;
    
    figure(5);
    hold on;
    plot(1:1:t, y(1,1:t));
    ylabel('State y_t');
    xlabel('Time (T)');
    pause(0.01);
    hold off;

    %%%%%%%%%%%%%%%%%%% MPC step end %%%%%%%%%%%%%%%%%%%

    u = v_mpc;

    % Measure state FIXME : t > 2
    if t > 2
        y(t+1) = b1*u(t) + b2*u(t-1) - a1*y(t) - a2*y(t-1) + w(t+1);
    end

    x(1,t+1) = y(t+1);

    if t >= 30 && mod(t, T_recalc) == 0
        %%%%%%%%%%%%%%%%%%% LSCR step start %%%%%%%%%%%%%%%%%%%
        result0 = LSCR_ARX((t-2)/r,N,M,y,u,w,D,0,r,1,theta0,theta1,2);
        result1 = LSCR_ARX((t-2)/r,N,M,y,u,w,D,1,r,1,theta0,theta1,2);

        result_common = result0.*result1;

        result3 = LSCR_ARX((t-2)/r,N,M,y,u,w,D,0,r,0,theta0,theta1,1);

        result_ab = zeros(N, N);

        % Claculate unknown parameters a1 and b1 from theta0 and theta1
        for i=1:N
            for j=1:N
                t0 = 2*i/N;
                t1 = b2 - (2*j/N) * t0;
                if (t1 > 0) && (t1 < 2)
                    k = round(N*t1/2);
                    if k > 0 && result_common(i,k) == 1
                        result_ab(i,j) = 1;
                    end
                end
            end
        end

        Omega_AB{3} = result_ab;

        figure(3);
        hold on;
        pcolor(theta, theta, result_common);
        xlabel('theta0 (= b0)');
        ylabel('theta1 (= b1 - a0 * b0)');
        plot(theta0, theta1, 'r*');
        pause(0.01);
        hold off;
        
        figure(4);
        hold on;
        pcolor(theta, theta, result_ab);
        xlabel('a0');
        ylabel('b0');
        pause(0.01);
        hold off;
        %%%%%%%%%%%%%%%%%%% LSCR step end %%%%%%%%%%%%%%%%%%%
    end
end