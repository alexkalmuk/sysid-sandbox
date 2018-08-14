function result = MPCCostFunc2(u,A_arr,B_arr,W_arr,x0,N,M)
    res(1:M) = cvx(zeros(1,M));

    for i=1:M % scenarios
        % Get i'th scenario
        A = A_arr(:,:,i);
        B = B_arr(:,:,i);
        W = W_arr(:,:,i);

        x = x0;

        % The main idea is to optimize along different paths with a common source,
        % which correspond to selected scenarios. THe similar
        % idea is called Policy Optimized MPC, and can be investigated
        % in more detail here -
        % "Min–Max Feedback Model Predictive Control for Constrained Linear
        % Systems" by P. O. M. Scokaert and D. Q. Mayne.

        % We use branches u(1), u(2), u(3),...
        % u(1), u(2+N-1),u(2+N),...
        % u(1), u(2+(N-1)*2), u(2+(N-1)*2+1), ...
        %
        z = A * x' + B' * u(1) + [W(1) W(1)]';
        x = z';
        res(i) = res(i) + x(1)^2 + u(1)^2;

        index = (N-1)*(i-1);

        for j=2:N % horizon
            % We use x(1)^2 instead of sum(x.^2) because we are interested
            % only in x(1), which is corresponds y = C * x, where C = [1 0
            % 0 0 ... 0].
            z = A * x' + B' * u(index+j) + [W(j) W(j)]';
            x = z';
            res(i) = res(i) + x(1)^2 + u(index+j)^2;
        end
        res(i) = res(i) + 100 * (x(1)^2);
    end
    result = max(res);
end