clear all
clc

N = 10^4;

e = 1*randn(N, 1);
e_u = randn(N, 1);

G_eu = tf([1 0.2], [1 -0.1 -0.12], -1);
G_ey = tf([1 -1 0.2], [1 -1.5 0.7], -1);
G_uy = tf([1 0.5], [1 -1.5 0.7], -1);

u = lsim(G_eu, e_u);
y = lsim(G_uy, u) + lsim(G_ey, e);

C_inv = tf([1], [1 -1 0.2], -1);
y_f = lsim(C_inv, y);
u_f = lsim(C_inv, u);

phi = zeros(N,4);
phi(1,:) = [0, 0, 0, 0];
phi(2,:) = [y(1), 0, u(1), 0];
for k = 3:50
    phi(k,:) = [-y_f(k-1), -y_f(k-2), u_f(k-1), u_f(k-2)];
end
theta = phi\y_f;
disp(theta);

G_est = tf([theta(3) theta(4)], [1 theta(1) theta(2)], -1);
y_prdt = lsim(G_est, u) + lsim(G_ey, e);
figure(1)
plot(y(1001:2000))
hold on
plot(y_prdt(1001:2000))
