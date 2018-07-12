function [ v ] = solve_rmpc_cop(A_scen,B_scen,W_scen,C_u,x_cur,N,S)
% Solves Convex Optimization problem for RMPC
% A_scen, B_scen,W_scen are the random selected scenarios
% S - the number of scenarios
% N - horizon length
% C_u - constraints

cvx_begin quiet
    variables c u(N)
    %expression x0(1,2);
    expression res(1,S);
    %x0 = x_cur_mcp;
    minimize ( c )
    subject to
        MPCCostFunc(u,A_scen,B_scen,W_scen,x_cur,N,S) <= c;
        for j=1:N
            C_u(1) <= u(j) <= C_u(2);
        end
cvx_end

% Return only the first u(1)
v = u(1);

end

