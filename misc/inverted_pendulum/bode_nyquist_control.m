%
% Based on 
% http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlFrequency
%

M = 0.5;
m = 0.2;
b = 0.1;
I = 0.006;
g = 9.8;
l = 0.3;
q = (M+m)*(I+m*l^2)-(m*l)^2;
s = tf('s');
P_pend = (m*l*s/q)/(s^3 + (b*(I + m*l^2))*s^2/q - ((M + m)*m*g*l)*s/q - b*m*g*l/q);

[zeros poles] = zpkdata(P_pend,'v')

% Derived from https://en.wikipedia.org/wiki/Argument_principle :
% Wind number (at the point -1) is N - P, where N and P denote number
% of zeros of poles respectively.
K = 10;
C = K*(s + 1)^2/(s);

figure('Name', 'Bode and Nyquist for P_pend and P_pend * C', ...
    'Position', [10 500 800 600], 'Color', 'w');
figure(1);
hold on;
    subplot(2,2,1);
    bodeplot(P_pend, 'r');

    subplot(2,2,2);
    bodeplot(P_pend * C, 'r');

    subplot(2,2,3);
    nyquist(P_pend);

    subplot(2,2,4);
    nyquist(P_pend * C);
hold off;

% Finally, tune the obtained controller
K = 35;
C = K*(s + 1)^2/(s);
T = feedback(P_pend,C);

% Calculate cart position
P_cart = (((I+m*l^2)/q)*s^2 - (m*g*l/q))/(s^4 + (b*(I + m*l^2))*s^3/q - ((M + m)*m*g*l)*s^2/q - b*m*g*l*s/q);
T2 = feedback(1,P_pend*C)*P_cart;
T2 = minreal(T2);

figure('Name', 'Pendulum and cart positions due to impulse', ...
    'Position', [750 500 600 600], 'Color', 'w');
figure(2);
hold on;
    subplot(2,1,1);
    t = 0:0.01:10;
    impulse(T,t), grid
    title({'Response of Pendulum Position to an Impulse Disturbance';'under Closed-loop Control'});

    subplot(2,1,2);
    t = 0:0.01:10;
    impulse(T2, t), grid
    title({'Response of Cart Position to an Impulse Disturbance';'under Closed-loop Control'});
hold off;