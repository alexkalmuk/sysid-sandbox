% Robust randomized MPC.
%   C - parameters of the initial system
%   C_u - C_u are constraints for the control u -- interval
%   Omega_AB - set of possible parameter values (uncertainties)
%   Omega_W - set of possible noise values
%   values
%   v_all - previos controls, which are, used ONLY for plotting
%   x_all - previos states, which are used ONLY for plotting
%   w - noise in the initial ARX system
%   x0 - initial state
%   N - MPC prediction horizon
%   M - scenarios count for Scenario Approach
%   t0 - start time instant
%
% Return: Array of contols - v, and the last state x_final.
function [v, x_final, res] = RMPC(cop_solver, C, C_u, C_y, Omega_AB,...
                            Omega_W, w, v_all, x_all, N, S, t0)
    A = [-C(1)   1;
         -C(2)   0];

    B = [C(3) C(4)];
    
    W = w' * [-C(1) -C(2)];
    
    % Selected scenarios for parameters in A and B
    % They are randomly chosen values from the set Omega_AB{3},
    % which is N x N epsilon set (N is Omega_AB{2}) based on
    % Omega_AB{1}.
    Omega_bounds = Omega_AB{1};
    Omega_N = Omega_AB{2};
    Omega_set = Omega_AB{3};
    
    found = 0;
    scen_A = zeros(1, S);
    scen_B = zeros(1, S);
    while found < S
        i = randi(Omega_N);
        j = randi(Omega_N);
        
        if Omega_set(i,j) == 1
           scen_a = Omega_bounds(1,1) + (Omega_bounds(1,2) - Omega_bounds(1,1))* i / Omega_N;
           scen_b = Omega_bounds(1,1) + (Omega_bounds(2,2) - Omega_bounds(2,1))* j / Omega_N;
           found = found + 1;
           scen_A(found) = scen_a;
           scen_B(found) = scen_b;
        end
    end

    % Selected scenarios for noise
    scen_W = zeros(1,N,S);
    for i=1:S
        scen_W(:,:,i) = Omega_W(1) + (Omega_W(2) - Omega_W(1)).*rand(1, N);
    end
    
    % Generate A and B under uncertainties
    A_mpc = zeros(2,2,S);
    B_mpc = zeros(1,2,S);
    W_mpc = scen_W;
    for i=1:S
        A_tmp = [-scen_A(i) A(1,2);
                 A(2,1)     A(2,2)];
        B_tmp = [scen_B(i)  B(2)];

        A_mpc(:,:,i) = A_tmp;
        B_mpc(:,:,i) = B_tmp;
    end

    x_cur = x_all(:,t0)';

    v = cop_solver(A_mpc,B_mpc,W_mpc,C_u,C_y,x_cur,N,S);

    if isnan(v)
        fprintf(['Error: Convex Optimization Problem cannot be solved!\n',...
                '    Try to increase y_t constraints C_y\n']);
        v = v_all;
        x_final = x_all;
        res = -1; % Error
        return;
    end

    % Calculate next x_cur=x_{t+1} using obtained u(1)
    z = A * x_cur' + B' * v + W(i, :)';
    x_cur = z';

    figure(2);
    hold on;

    subplot(2,1,1);
    plot(1:(t0), horzcat(v_all, v), 'Color', 'blue');
    xlabel('Time (T)');
    ylabel('Control (v_t)');

    subplot(2,1,2);
    plot(1:t0, x_all(1,1:t0), 'Color', 'blue');
    xlabel('Time (T)');
    ylabel('Approximated (x_t(1))');
    hold off;
    
    v = horzcat(v_all, v);
    x_final = horzcat(x_all, x_cur');
    res = 0; % Success
end