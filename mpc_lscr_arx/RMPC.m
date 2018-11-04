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

    [A_mpc, B_mpc, W_mpc] = gen_scenarios(A,B,Omega_AB,Omega_W,N,S);

    x_cur = x_all(:,t0)';

    v = cop_solver(A_mpc,B_mpc,W_mpc,C_u,C_y,x_cur,N,S);

    if isnan(v)
        fprintf(['Error: Convex Optimization Problem cannot be solved!\n',...
                '    Try to increase y_t State Constraints\n']);
        disp('x_cur:');
        disp(x_cur);

        v = v_all;
        x_final = x_all;
        res = -1; % Error
        return;
    end

    % Calculate next x_cur=x_{t+1} using obtained u(1)
    z = A * x_cur' + B' * v + W(1, :)';
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