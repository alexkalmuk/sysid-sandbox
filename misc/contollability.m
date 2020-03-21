% ARX system
a1 = 1.4;
a2 = 1;
b1 = 0.6;
b2 = 1.6;

A = [-a1   1; 
     -a2   0];
 
B = [b1 b2]';
C = [1 0];

Ctrl = ctrb(A,B);

uncontrollable = length(A) - rank(Ctrl);

if uncontrollable == 0
    fprintf('The system is controllable\n');
else
    fprintf('The system is uncontrollable\n');
end

sys_ss = ss(A,B,C,0);
tf_ss = tf(sys_ss);

fprintf('State space system:\n');
sys_ss

fprintf('Transfer function:\n');
tf_ss

fprintf('Check if system is stable (Bode plot make sence only if the system is stable):');
isstable(sys_ss)

sys_ss_from_tf = ss(tf_ss);

figure('Name', 'Bode plots', 'Position', [10 500 600 400], 'Color', 'w');
figure(1);
bodeplot(tf_ss, 'r', sys_ss_from_tf, 'y--');

figure('Name', 'Nyquist plot', 'Position', [800 500 600 400], 'Color', 'w');
figure(2);
nyquist(sys_ss);

% Just feed sin to the system and see what happen
t = 0:0.1:20;
u = sin(0.5*t);
figure('Name', 'Simulating tf response', 'Position', [650 10 600 400], 'Color', 'w');
figure(3);
lsim(tf_ss,u,t);