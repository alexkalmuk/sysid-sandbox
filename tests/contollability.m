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
ts_ss = tf(sys_ss);

bodeplot(sys_ss)