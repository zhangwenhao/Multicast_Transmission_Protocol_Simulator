clc;
clear;
close all;

% Target rate
R = 0.1;
% Channel usage count
N = 100;

total_size=12;
looptimes= 10;
 
Q_VALUE = zeros(1, total_size, 'double'); 
APP_VALUE = zeros(1, total_size, 'double'); 
num_time = zeros(1, total_size, 'double');

varsigma = sqrt(N / (2 * pi * (exp(2 * R) - 1)));
vartheta = exp(R) - 1;

temp_sub=vartheta-1/(2*varsigma);
temp_add=vartheta+1/(2*varsigma); 
temp_sub=0.05;
temp_add=0.12;

q_func=@(x) 0.5 * erfc(x / sqrt(2));%Q function 与erfc 函数关系Q=1/2erfc(x/sqrt(2));

Q_VALUE_FUN = @(gamma_ag) q_func( (sqrt(N)*(1+gamma_ag)*(log(1+gamma_ag)-R))/(sqrt(gamma_ag*(gamma_ag+2))));

APP_VALUE_FUN = @(gamma_ag_temp )-varsigma*(gamma_ag_temp-vartheta)+0.5;

for index=0:total_size
    gamma_ag = temp_sub + (temp_add-temp_sub)/total_size * index;
    Q_VALUE(index+1) = Q_VALUE_FUN(gamma_ag);
    APP_VALUE(index+1) =  APP_VALUE_FUN(gamma_ag);
    num_time(index+1)=  gamma_ag;
end

figure(1); 
plot(num_time, Q_VALUE,'Color','k','LineStyle','-','Marker','square','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(num_time),'linewidth',1);
hold on;
plot(num_time, APP_VALUE,'Color','k','LineStyle','--','Marker','o','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(num_time),'linewidth',1);
 
xlabel('SNR between UAV and k-th user, $\gamma_{ak}$' ,'Interpreter','latex');
ylabel('Decoding error probability','Interpreter','latex');
legend('Q-Function','Approximation','Location','best');
box on;
grid on;