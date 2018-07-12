function result = MPCCostFunc(u,A_arr,B_arr,W_arr,x0,N,M)
    res(1:M) = cvx(zeros(1,M));

    for i=1:M % scenarios
        % Get i'th scenario
        %A = A_arr(:,:,i);
        %B = B_arr(:,:,i);
        % FIXME Remove that hardcoded values and change back to the
        % presented above
        A = [-1.4 1;
             -1   0];
        B = [0.6 1.6];

        W = W_arr(:,:,i);

        x = x0;

        for j=1:N % horison
            % We use x(1)^2 instead of sum(x.^2) because we are interested
            % only in x(1), which is corresponds y = C * x, where C = [1 0
            % 0 0 ... 0].
            z = A * x' + B' * u(j) + [W(j) W(j)]';
            x = z';
            res(i) = res(i) + x(1)^2 + u(j)^2;
        end
        res(i) = res(i) + 100 * (x(1)^2);
    end
    result = max(res);
end