function result = LSCR_ARX(T, N, M, y, u, w, D, s, r, l, theta0, theta1, dim, bounds)
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
                    y_estim(t,i,j) = u(t+1) * i_val;
                    z_estim = D(t+1) * i_val + D(t);
                else
                    y_estim(t,i,j) = u(t+1) * i_val + u(t) * j_val;
                    z_estim = D(t+1) * i_val + D(t) * j_val;
                end
                e(t,i,j) = (y(r*t+l) - y_estim(t,i,j)) - (getZ(D(t+1), D(t), w(t)) - z_estim);
                f(t,i,j) = sign((D(t+1-s) + u(t+1-s)) * e(t,i,j));
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

    function d = getZ(D1, D2, w)
        d = theta0 * D1 + theta1 * D2 + w;
    end
end