clc;
clear;
close all;

% Set the random seed to ensure consistent random number generation each time
rng(5);

% Cluster data (number of points, Poisson cluster count, clustering probability, radius)
radius = 1000;
numberNodes = 50;
centersSize = 5;
clusterParameter = 0.3;

[data, k] = poisson_cluster(numberNodes, centersSize, clusterParameter, radius); 

loop = true;

warden = [0, 75];
%warden = [0, -50];
%warden = [0, -70];

% Plotting the initial state
%plot_gu(warden, data, radius);

% Radius for guard zone
r_w = 150;

% Perform K-means++ clustering
[idx, ctr, wdx] = group_k_means(data, warden, r_w);

% Plot the results of K-means++ clustering
plot_mg(warden, wdx, ctr, radius, r_w, idx);

% Perform uniform radius grouping
[MBSLocations, finalRadius, sortedWdx] = group_uniform_radius(data, warden, r_w, radius/2);

% Plot the results of uniform radius grouping
plot_uniform_group(warden, data, MBSLocations, sortedWdx, r_w, radius, finalRadius);

% Number of groups for K-means++ clustering
k_kmeans = size(ctr, 1);
% Number of GUs within guard zone for K-means++
g_kmeans = size(wdx, 1);

% Number of groups for MBS
k_mbs = size(MBSLocations, 1);
% Number of GUs within guard zone for MBS
g_mbs = size(sortedWdx, 1);

% Target rate
R = 0.1;
% Channel usage count
N = 100;
% Data size
M = 1 * 1024 * 1024;

% Noise variance
sigma_k = fun_db_to_math(-70);
sigma_w = sigma_k;
% UAV height
h = 500;
% Channel coefficient
beta = fun_db_to_math(-40);
% Transmission probability
p1 = 0.5; 

TOTAL_TIME_KMEANS = zeros(1, k_kmeans + g_kmeans, 'double');
num_time_KMEANS = zeros(1, k_kmeans + g_kmeans, 'double');

TOTAL_TIME_MBS = zeros(1, k_mbs + g_mbs, 'double');
num_time_MBS = zeros(1, k_mbs + g_mbs, 'double');

% covertness constraint
epslon = 0.1;
varsigma = sqrt(N / (2 * pi * (exp(2 * R) - 1)));
vartheta = exp(R) - 1;
eta_k = @(gamma_ak_temp) -varsigma * (gamma_ak_temp - vartheta) + 0.5;
C_k = @(eta_k_temp) N * R * (1 - eta_k_temp);

% Calculate transmission time for GUs within guard zone - K-means++
for i = 1:g_kmeans
    h_ak = beta / power(h, 2);
    h_aw = beta / (power(h, 2) + power(norm(wdx(i, 1:2) - warden), 2)); 
    p_ak = 4 * epslon * sigma_w * sqrt(2 / N) / h_aw;
    gamma_ak = p_ak * h_ak / sigma_k;
    R_eta_k = eta_k(gamma_ak);
    R_C_k = C_k(R_eta_k);
    TOTAL_TIME_KMEANS(i) = M / (p1 * R_C_k);
    num_time_KMEANS(i) = i;
end

% Calculate transmission time for K-means++ clusters
for j = 1:k_kmeans
    h_ak = beta / (power(h, 2) + power(ctr(j, 3), 2));
    h_aw = beta / (power(h, 2) + power(norm(ctr(j, 1:2) - warden), 2)); 
    p_ak = 4 * epslon * sigma_w * sqrt(2 / N) / h_aw;
    gamma_ak = p_ak * h_ak / sigma_k;
    R_eta_k = eta_k(gamma_ak);
    R_C_k = C_k(R_eta_k);
    TOTAL_TIME_KMEANS(j + g_kmeans) = M / (p1 * R_C_k);
    num_time_KMEANS(j + g_kmeans) = j + g_kmeans;
end

% Calculate transmission time for GUs within guard zone - MBS
for i = 1:g_mbs
    h_ak = beta / power(h, 2);
    h_aw = beta / (power(h, 2) + power(norm(sortedWdx(i, 1:2) - warden), 2)); 
    p_ak = 4 * epslon * sigma_w * sqrt(2 / N) / h_aw;
    gamma_ak = p_ak * h_ak / sigma_k;
    R_eta_k = eta_k(gamma_ak);
    R_C_k = C_k(R_eta_k);
    TOTAL_TIME_MBS(i) = M / (p1 * R_C_k);
    num_time_MBS(i) = i;
end

% Calculate transmission time for MBS clusters
for j = 1:k_mbs
    h_ak = beta / (power(h, 2) + power(MBSLocations(j, 3), 2));
    h_aw = beta / (power(h, 2) + power(norm(MBSLocations(j, 1:2) - warden), 2)); 
    p_ak = 4 * epslon * sigma_w * sqrt(2 / N) / h_aw;
    gamma_ak = p_ak * h_ak / sigma_k;
    R_eta_k = eta_k(gamma_ak);
    R_C_k = C_k(R_eta_k);
    TOTAL_TIME_MBS(j + g_mbs) = M / (p1 * R_C_k);
    num_time_MBS(j + g_mbs) = j + g_mbs;
end

% Plot results for K-means++ clustering
figure(1); 
barWidth = 0.5;
b1=bar(num_time_KMEANS, TOTAL_TIME_KMEANS, barWidth, 'FaceColor', 'k', 'EdgeColor', 'k'); 
xlabel('Serial number of MG');
ylabel('Transmission time (ms)');
xticks(1:1:g_kmeans + k_kmeans);
legend('theoretical results', 'Location', 'northeast');
box on;
grid on;
hold off;

% Plot results for MBS
figure(2); 
barWidth = 0.5;
b3=bar(num_time_MBS, TOTAL_TIME_MBS, barWidth, 'FaceColor', 'k', 'EdgeColor', 'k'); 
xlabel('Serial number of MG');
ylabel('Transmission time (ms)');
xticks(1:1:g_mbs + k_mbs);
legend('theoretical results', 'Location', 'northeast');
box on;
grid on;
hold off;
