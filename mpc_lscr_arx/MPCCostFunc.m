function result = MPCCostFunc(u,A_arr,B_arr,W_arr,x0,N,M)
    res(1:M) = cvx(zeros(1,M));
    x = x0;
    for i=1:M % scenarios
        % Get i'th scenario
        A = A_arr(:,:,i);
        B = B_arr(:,:,i);
        W = W_arr(:,:,i);
        for j=1:N % horison
            % We use x(1)^2 instead of sum(x.^2) because we are interested
            % only in x(1), which is corresponds y = C * x, where C = [1 0
            % 0 0 ... 0].
            res(i) = res(i) + 4 * x(1)^2 + u(j)^2;
            z = A * x' + B' * u(j) + W(j);
            x = z';
        end
    end
    result = max(res);
end