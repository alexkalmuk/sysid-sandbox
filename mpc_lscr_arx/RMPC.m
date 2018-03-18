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
%   T - final time instant 
%
% Return: Array of contols - v, and the last state x_final.
function [v, x_final] = RMPC(C, C_u, Omega_AB, Omega_W, w, v_all, x_all, N, S, t0, T_final)
    A = [-C(1)   1;
         -C(2)   0];

    B = [C(3) C(4)];
    
    W = w' * [-C(1) -C(2)];
    
    % The time interval from the start to the end
    T = T_final - t0 + 1;
    
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
    
    disp(scen_A);
    disp(scen_B);

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
    
    v = zeros(1,T);
    x = zeros(2,T+1);
    x_cur = x_all(:,t0)';

    for i=1:T
        x(:,i) = x_cur';
        
        % Calculate current u_t = u(1) under horizon N_mpc
        cvx_begin quiet
            variables c u(N)
            %expression x0(1,2);
            expression res(1,S);
            %x0 = x_cur_mcp;
            minimize ( c )
            subject to
                MPCCostFunc(u,A_mpc,B_mpc,W_mpc,x_cur,N,S) <= c;
                for j=1:N
                    C_u(1) <= u(j) <= C_u(2);
                end
        cvx_end

        % Save obtained u(1)
        v(i) = u(1);

        % Calculate next x_cur=x_{t+1} using obtained u(1)
        z = A * x_cur' + B' * u(1) + W(i, :)';
        x_cur = z';
        
        figure(2);
        
        subplot(2,1,1);
        plot(1:(t0+i-1), horzcat(v_all, v(1:i)));
        xlabel('Time (T)');
        ylabel('Control (v_t)');
        
        subplot(2,1,2);
        plot(1:(t0+i), horzcat(x_all(1,1:t0), x(1,1:i)));
        xlabel('Time (T)');
        ylabel('State (y_t)');
        
        pause(0.01)
    end
    
    hold off;
    
    v = horzcat(v_all, v);
    x_final = horzcat(x_all, x);
end