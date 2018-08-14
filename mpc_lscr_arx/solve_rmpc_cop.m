function [ v ] = solve_rmpc_cop(A_scen,B_scen,W_scen,C_u,C_y,x_cur,N,S)
% Solves Convex Optimization problem for RMPC
% A_scen, B_scen,W_scen are the random selected scenarios
% S - the number of scenarios
% N - horizon length
% C_u, C_y - constraints

cvx_begin quiet
    variables c u(N)
    expression res(1,S);
    minimize ( c )
    subject to
        MPCCostFunc(u,A_scen,B_scen,W_scen,x_cur,N,S) <= c;
        % Control constraints
        for j=1:N
            C_u(1) <= u(j) <= C_u(2);
        end
        % State constraints
        for j=1:S
            A = A_scen(:,:,j);
            B = B_scen(:,:,j);
            W = W_scen(:,:,j);
            x = x_cur';
            for i=1:N
                for k = 1:N
                    x = cvx(A * x) + B' * u(i) + [W(k) W(k)]';
                    C_y(1) <= x(1) <= C_y(2);
                end
            end
        end
cvx_end

% Return only the first u(1)
v = u(1);

end

