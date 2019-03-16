% Based on http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlStateSpace

M = 0.5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;

p = I*(M+m)+M*m*l^2; %denominator for the A and B matrices

A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];
B = [     0;
     (I+m*l^2)/p;
          0;
        m*l/p];
C = [1 0 0 0;
     0 0 1 0];
D = [0;
     0];

% Assume that our reference signal is [r, 0] = [0.2, 0],
% that is we want to control cart position to 0.2 meters rigth
% adn angle=0.
r = 0.2;

sys_ss = ss(A,B,C,D);

ob = obsv(sys_ss);
observability = rank(ob);

ob = obsv(sys_ss);
unobservable = length(A) - rank(ob);

fprintf('State space system of inverted pendulum:\n');
sys_ss

if unobservable == 0
    fprintf('The system is observable\n');
else
    fprintf('The system is ununobservable\n');
end

% Use LQR
Q = C'*C;
Q(1,1) = 5000; % tunable
Q(3,3) = 100;  % tunable
R = 1;
K = lqr(A,B,Q,R);

% Pole placement design
% So we place poles at positions we want.
% Just to choose some negative for stability reasons.
P = [-40 -41 -42 -43];
L = place(A', C', P)';

% Assume that our reference signal is [r, 0] = [0.2, 0].
%
% Adding precompensation. Look at Astrem book 6.2 Stabilization by State Feedback
% "State Space Controller Structure"
%
% x_e = -(A - B*K)*B*k_r*r;
% y_e = C*x_e;
%
% y_e = [r, 0] = -C(A - B*K)*B*k_r*r;
%
M = -C*inv(A - B*K)*B;

% I means, [r, 0] = M*k_r*r. So, M(1)*k_r = 1 => k_r = 1/M(1)
k_r = 1 / M(1);

Ac = [(A-B*K) (B*K);
       zeros(length(A)) (A-L*C)];
Bc = [B*k_r;
       zeros(size(B))];
Cc = [C zeros(size(C))];
Dc = [0;0];

sys_est_cl = ss(Ac,Bc,Cc,Dc);

t = 0:0.01:5;
r_sig = r*ones(size(t));
[y,t,x]=lsim(sys_est_cl,r_sig,t);
[AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
set(get(AX(1),'Ylabel'),'String','Cart position (m)')
set(get(AX(2),'Ylabel'),'String','Pendulum angle (radians)')
title('Step Response with Observer-Based State-Feedback Control')