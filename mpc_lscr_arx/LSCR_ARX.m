function result = LSCR_ARX(T, N, M, y, u, w, D, s, theta0, theta1, dim, bounds)
    a = bounds(1,1);
    b = bounds(1,2);
    c = bounds(2,1);
    d = bounds(2,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT
    result = zeros(N,N);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    y_estim = zeros(T,N,N);
    e = zeros(T,N,N);
    f = zeros(T,N,N);
    g = zeros(N,N);
    h = zeros(M,T);
    for i = 1:M
       h(i,1:T) = (rand(1, T) >= 0.5);
    end

    
    for t = 1:T
        for i = 1:N
            for j = 1:N
                i_val = a + (b-a)*i/N;
                j_val = c + (d-c)*j/N;

                if dim == 1
                    y_estim(t,i,j) = u(2*t-1) * i_val;
                    z_estim = D(2*t-1) * i_val;
                    z = D(2*t-1) * theta1 + w(2*t);
                else
                    y_estim(t,i,j) = u(2*t-1) * j_val + u(2*t) * i_val;
                    z_estim = D(2*t-1) * j_val + D(2*t) * i_val;
                    z = D(2*t-1) * theta0 + D(2*t) * theta1;
                end
                e(t,i,j) = (y(2*t+dim-1) - y_estim(t,i,j)) - (z - z_estim);
                f(t,i,j) = sign((D(2*t-s) + u(2*t-s)) * e(t,i,j));
            end
        end
    end

    for i = 1:M
        z = zeros(N,N);
        for j = 1:N
            for k = 1:N
                z(j,k) = h(i, 1:T) * f(1:T,j,k);
            end
        end
        g(i,1:N,1:N) = z;
    end

    for j = 1:N
        for k = 1:N
            l=0;
            r=0;
            for i=1:M
                if g(i,j,k) > 0
                    l=l+1;
                elseif g(i,j,k) < 0
                    r=r+1;
                end
                if l>2 && r>2
                    result(j,k) = 1;
                end
            end
        end
    end
end