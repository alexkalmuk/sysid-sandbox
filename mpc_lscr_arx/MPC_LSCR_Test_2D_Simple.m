%
% MPC + LSCR for ARX models
%

%%%%% System setup. a1 and b1 are unknown. %%%%%
a1 = 1.4;
a2 = 1;
b1 = 0.6;
b2 = 1.6;
% Initial state
x0 = ones(1,2) * 7;
% Time
T = 100;
% Time instants of kind k * T_recalc when a new set is calculated using
% LCSR_ARX
T_recalc = 20;

% Epsilon-set dimenstions. That is if the set of pissible parameters is
% rectangle with size a, b than precision is a / N, b / N.
N = 50;
% Omega is a set of the model uncertainties.
% a1 in [Omega_AB{1}(1,1), Omega(1,2)]
% b1 in [Omega(2,1), Omega(2,2)]
% w in [Omega(3,1), Omega(3,2)]
% But initially it is equal to the whole Omega.
Omega_AB = {[0 2; 0 2], N, ones(N, N)};
Omega_W = [0 10];
% C_u are constraints for the control u -- interval
C_u = [-1 5];

% Noise is unknown but bounded and located within [Omega(3,1), Omega(3,2)]
mu = 4;
sigma = 7;
w = normrnd(mu,sigma,[1 T]);
w(w > Omega(3,2)) = mu;
w(w < Omega(3,1)) = mu;

%%%%% MPC setup %%%%%
N_mpc = 3; % horizon
M_mpc = 10; % number of scenarios

%%%%% LSCR_ARX Setup %%%%%
M = 100;
r = 2;
theta = 0:(2/N):(2-2/N);

%%%%% Below we do the main part - MPC + LSCR_ARX %%%%%
C = [a1 a2 b1 b2];

iterations = round(T / T_recalc);
T0 = 1;
x = zeros(2, 1);
x(:, 1) = x0';

figure('Name', 'RMPC', 'Position', [10 500 600 290]);
figure(2);
figure('Name', 'Controls u_t and u_t + D', 'Position', [10 10 600 290]);
figure(3);
figure('Name', 'Result in terms of theta0 and theta1', 'Position', [650 500 600 290]);
figure(4);
figure('Name', 'Final result for the current step', 'Position', [650 10 600 290]);
figure(5);

v_mpc = [];
for iter=1:iterations
    T = iter * T_recalc;
    
    %%%%%%%%%%%%%%%%%%% MPC step start %%%%%%%%%%%%%%%%%%%

    [v_mpc, x] = RMPC(C, C_u, Omega_AB, Omega_W, w, v_mpc, x, N_mpc, M_mpc, T0, T);

    %%%%%%%%%%%%%%%%%%% MPC step end %%%%%%%%%%%%%%%%%%%
    
    D = sign(randn(1,T));
    D(D == 0) = 1;
    D = 25 * D;
    
    figure(3);

    plot(1:1:T, v_mpc); hold on;
    plot(1:1:T, v_mpc + D); hold on;
    plot(1:1:T, zeros(1,T));
    ylabel('v_t (blue), \Delta_t (red)');
    xlabel('Time (T)');
    hold off;

    %%%%%%%%%%%%%%%%%%% LSCR step start %%%%%%%%%%%%%%%%%%%
    theta0 = b1;
    theta1 = (b2 - a1 * b1);

    %u = normrnd(20,3,[1 3*T]);
    u = v_mpc;

    y = zeros(T,N,N);
    for t=3:T-1
        y(t,1:N,1:N) = (b1*u(t+1) + b2*u(t) - a1*y(t-1,1,1) - a2*y(t-2,1,1) + w(t)) * ones(N,N);
    end

    result0 = LSCR_ARX((T-2)/r,N,M,y,u,w,D,0,r,1,theta0,theta1,2);
    result1 = LSCR_ARX((T-2)/r,N,M,y,u,w,D,1,r,1,theta0,theta1,2);

    result_common = result0.*result1;

    result3 = LSCR_ARX((T-2)/r,N,M,y,u,w,D,0,r,0,theta0,theta1,1);

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

    %pcolor(theta, theta, result0); figure;
    %pcolor(theta, theta, result1); figure;
    figure(4);
    pcolor(theta, theta, result_common);
    xlabel('theta0 (= b1)');
    ylabel('theta1 (= b2 - a1 * b1)');
    figure(5);
    pcolor(theta, theta, result_ab);
    xlabel('a1');
    ylabel('b1');
    %pcolor(theta, theta, result3);
    %%%%%%%%%%%%%%%%%%% LSCR step end %%%%%%%%%%%%%%%%%%%
    
    T0 = T0 + T_recalc;
end
