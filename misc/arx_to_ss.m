a1 = 1.4;
a2 = 1;
b1 = 0.6;
b2 = 1.6;

A = [-a1   1; 
     -a2   0];
 
B = [b1 b2];
C = [1 0];

T = 20;
N = 2;

x = [1 1];

u = ones(1,T+5);
y = zeros(1,T+5);
out = zeros(1,T+5);

for i=1:5
    z = A * x' + B' * u(i);
    y(i) = x(1);
    x = z';
end

for i=5:T
    z = A * x' + B' * u(i);
    out(i+1) = x(1);
    x = z';
    y(i+1) = -a1*y(i) - a2*y(i-1) + b1*u(i) + b2*u(i-1);
end

fprintf('Result:\n');
fprintf('ARX:\n');
disp(y(10:18));
fprintf('State space:\n');
disp(out(10:18));