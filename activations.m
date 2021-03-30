close all;
clear;
clc;
x = -5:0.01:5;
%step function
step = x>0;
%sigmoid function
sigmoid = 1./(1+exp(-x));
%relu fucntion
relu = x+0.05;
relu(relu<0) = 0;
%leaky relu function
leaky_relu = x;
leaky_relu(x<0) = 0.1*x(x<0);
%tanh function
tan_h = 2./(1+exp(-2*x)) -1;

figure(1);
plot(x, step, 'linewidth',2);
hold on;
plot(x, sigmoid, 'linewidth',2);
hold on;
plot(x, relu, 'linewidth',3);
hold on;
plot(x, leaky_relu, 'linewidth',2);
hold on;
plot(x, tan_h, 'linewidth',2);

ylim([-1.5,2]);
xlim([-5,5]);
ax = gca;
ax.FontSize = 16; 
xlabel('X (input)','FontSize',20)
ylabel('f(X) (output)','FontSize',20)
legend({"heaviside step", "sigmoid","relu", "leaky relu, a = 0.01", "tanh"}, 'fontsize', 22);


x = -5:0.01:5;
%step function derivative
step = (x==0).*2;
%sigmoid function derivative (f(x)*(1-f(x))
sigmoid = (1./(1+exp(-x))).*(1-1./(1+exp(-x)));
%relu fucntion derivative
relu = x-0.05;
relu = relu>=0;
%leaky relu function derivative
leaky_relu = x;
leaky_relu(x>=0) = 1;
leaky_relu(x<0) = 0.1;
%tanh function derivative (1-f(x)^2)
tan_h = 1- (2./(1+exp(-2*x)) -1).^2;

figure(2);
plot(x, step, 'linewidth',2);
hold on;
plot(x, sigmoid, 'linewidth',2);
hold on;
plot(x, relu, 'linewidth',3);
hold on;
plot(x, leaky_relu, 'linewidth',2);
hold on;
plot(x, tan_h, 'linewidth',2);

ylim([-0.5,3]);
xlim([-5,5]);
ax = gca;
ax.FontSize = 16; 
xlabel('X (input)','FontSize',20)
ylabel('f(X) (output)','FontSize',20)
legend({"heaviside step", "sigmoid","relu", "leaky relu, a = 0.01", "tanh"}, 'fontsize', 22);

