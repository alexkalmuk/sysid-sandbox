function [ A_out, B_out, W_out ] = gen_scenarios(A, B, Omega_AB, Omega_W, N, S)   
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
    A_out = zeros(2,2,S);
    B_out = zeros(1,2,S);
    W_out = scen_W;
    for i=1:S
        A_tmp = [-scen_A(i) A(1,2);
                 A(2,1)     A(2,2)];
        B_tmp = [scen_B(i)  B(2)];

        A_out(:,:,i) = A_tmp;
        B_out(:,:,i) = B_tmp;
    end
end

