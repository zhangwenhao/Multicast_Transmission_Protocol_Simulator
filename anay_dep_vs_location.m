clc;
clear;
close all;

parpool(5);

% Radius of the area
radius = 1000;

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
% Covertness constraint
epsilon_001 = 0.2;
epsilon= 0.1;

lo_uncertainty=10;

looptimes=9;

% DEP
MG_DEP_001=zeros(1, looptimes+1, 'double');
UNIFORM_DEP_001=zeros(1, looptimes+1, 'double');
MG_DEP=zeros(1, looptimes+1, 'double');
UNIFORM_DEP=zeros(1, looptimes+1, 'double');

% Transmission time
MG_TIME_001=zeros(1, looptimes+1, 'double');
UNIFORM_TIME_001=zeros(1, looptimes+1, 'double');
MG_TIME=zeros(1, looptimes+1, 'double');
UNIFORM_TIME=zeros(1, looptimes+1, 'double');

% Total transmission time
MG_TIME_001_TOTAL=zeros(1, looptimes+1, 'double');
UNIFORM_TIME_001_TOTAL=zeros(1, looptimes+1, 'double');
MG_TIME_TOTAL=zeros(1, looptimes+1, 'double');
UNIFORM_TIME_TOTAL=zeros(1, looptimes+1, 'double');


loop_times = zeros(1,looptimes+1,'double');

dep = @(phi)0.5-0.5*sqrt(phi*N/4);

phi = @(p_ak,h_aw) p_ak*h_aw/sigma_w-log(1+p_ak*h_aw/sigma_w);

varsigma = sqrt(N/(2*pi*(exp(2*R)-1)));
vartheta = exp(R)-1;
eta_k = @(gamma_ak_temp) -varsigma*(gamma_ak_temp-vartheta)+0.5;
C_k = @(eta_k_temp) N*R*(1-eta_k_temp);

numberNodes=200;
centerszie=10;
cluster_parameter=0.5;
warden=[0,0];
warden_act=[0,0];
r_w = 150;

distributionLoop=20000;
groupLoop=8;

%parfor 
parfor index=0:looptimes
    lo_uncertainty=10*index;
    loop_times(index+1)=lo_uncertainty;
    distribution_dep_kmeans=0;
    distribution_dep_mbs=0;
    distribution_dep_kmeans_001=0;
    distribution_dep_mbs_001=0;
    distribution_Time_kmeans=0;
    distribution_Time_mbs=0;
    distribution_Time_kmeans_001=0;
    distribution_Time_mbs_001=0;
    distribution_Time_kmeans_total=0;
    distribution_Time_mbs_total=0;
    distribution_Time_kmeans_001_total=0;
    distribution_Time_mbs_001_total=0;
    warden=[0,lo_uncertainty];
    for d_l=1:distributionLoop
        fprintf('looptimes: %d, distributionloop: %d\n', index, d_l);
        centerszie_p = poissrnd(centerszie);
        [data, ~]=poisson_cluster(numberNodes, centerszie_p, cluster_parameter, radius); 
        %close all;
        %plot_gu(warden,data,radius);
        MG_CONS_dep_kmeans=0;
        MG_CONS_dep_mbs=0;  
        MG_CONS_dep_kmeans_001=0;
        MG_CONS_dep_mbs_001=0; 
        MG_CONS_Time_kmeans_total=0;
        MG_CONS_Time_mbs_total=0;  
        MG_CONS_Time_kmeans_001_total=0;
        MG_CONS_Time_mbs_001_total=0;
        MG_CONS_Time_kmeans=0;
        MG_CONS_Time_mbs=0;  
        MG_CONS_Time_kmeans_001=0;
        MG_CONS_Time_mbs_001=0; 
        for g_l=1:groupLoop
            if g_l>=groupLoop*3/4
                warden=[0,lo_uncertainty];
            elseif g_l>=groupLoop*2/4
                warden=[0,-lo_uncertainty];
            elseif g_l>=groupLoop/4
                warden=[lo_uncertainty,0];
            else
                warden=[-lo_uncertainty,0];
            end
            %fprintf('groupLoop: %d\n', g_l);
            [MBSLocations,finalRadius, sortedWdx]=group_uniform_radius(data,warden,r_w, radius/4);
            %plot_uniform_group(warden, data, MBSLocations, sortedWdx, r_w, radius, finalRadius);

            [idx, ctr, wdx] = group_k_means( data, warden, r_w, ceil(centerszie/2));
            %close all;
            %plot_mg(warden, wdx, ctr, radius, r_w, idx);

            MG_kmeans = [ctr; wdx];
            MG_mbs = [MBSLocations; sortedWdx];
    
            k_kmeans = size(MG_kmeans, 1);
            k_mbs = size(MG_mbs, 1);
            current_MG_dep_kmeans=0;
            current_MG_dep_mbs=0;
            current_MG_dep_kmeans_001=0;
            current_MG_dep_mbs_001=0;
            current_MG_Time_kmeans=0;
            current_MG_Time_mbs=0;
            current_MG_Time_kmeans_001=0;
            current_MG_Time_mbs_001=0;
            current_MG_Time_kmeans_total=0;
            current_MG_Time_mbs_total=0;
            current_MG_Time_kmeans_001_total=0;
            current_MG_Time_mbs_001_total=0;
            for i= 1:k_kmeans
                h_ak=beta/(power(h,2)+power(MG_kmeans(i,3),2));
                h_aw=beta/(power(h,2)+power(norm(MG_kmeans(i, 1:2) - warden_act),2)); 
                p_min=(vartheta-1/(2*varsigma))*sigma_k/h_ak;
                p_max=(vartheta+1/(2*varsigma))*sigma_k/h_ak;
                p_ak=4*epsilon*sigma_w*sqrt(2/N)/h_aw;
                p_ak_001=4*epsilon_001*sigma_w*sqrt(2/N)/h_aw;
                p_ak=max(p_ak,p_min);
                p_ak_001=max(p_ak_001,p_min);
                p_ak=min(p_ak,p_max);
                p_ak_001=min(p_ak_001,p_max);
                current_dep_kmeans=dep(phi(p_ak,h_aw));
                current_dep_kmeans_001=dep(phi(p_ak_001,h_aw));
                current_MG_dep_kmeans=current_MG_dep_kmeans + current_dep_kmeans;
                current_MG_dep_kmeans_001=current_MG_dep_kmeans_001 + current_dep_kmeans_001;
                current_time_kmeans=M/(p1*C_k(eta_k(p_ak*h_ak/sigma_k)));
                current_time_kmeans_001=M/(p1*C_k(eta_k(p_ak_001*h_ak/sigma_k)));
                current_MG_Time_kmeans=current_MG_Time_kmeans+current_time_kmeans;
                current_MG_Time_kmeans_001=current_MG_Time_kmeans_001+current_time_kmeans_001;
            end
            for i = 1:k_mbs
                h_ak=beta/(power(h,2)+power(MG_mbs(i,3),2));
                h_aw=beta/(power(h,2)+power(norm(MG_mbs(i, 1:2) - warden_act),2)); 
                p_min=(vartheta-1/(2*varsigma))*sigma_k/h_ak;
                p_max=(vartheta+1/(2*varsigma))*sigma_k/h_ak;
                p_ak=4*epsilon*sigma_w*sqrt(2/N)/h_aw;
                p_ak_001=4*epsilon_001*sigma_w*sqrt(2/N)/h_aw;
                p_ak=max(p_ak,p_min);
                p_ak_001=max(p_ak_001,p_min);
                p_ak=min(p_ak,p_max);
                p_ak_001=min(p_ak_001,p_max);
                current_dep_mbs=dep(phi(p_ak,h_aw));
                current_dep_mbs_001=dep(phi(p_ak_001,h_aw));
                current_MG_dep_mbs=current_MG_dep_mbs + current_dep_mbs;
                current_MG_dep_mbs_001=current_MG_dep_mbs_001 + current_dep_mbs_001;
                current_time_mbs=M/(p1*C_k(eta_k(p_ak*h_ak/sigma_k)));
                current_time_mbs_001=M/(p1*C_k(eta_k(p_ak_001*h_ak/sigma_k)));
                current_MG_Time_mbs=current_MG_Time_mbs+current_time_mbs;
                current_MG_Time_mbs_001=current_MG_Time_mbs_001+current_time_mbs_001;
            end
            MG_CONS_dep_kmeans=MG_CONS_dep_kmeans+current_MG_dep_kmeans/k_kmeans;
            MG_CONS_dep_mbs=MG_CONS_dep_mbs+current_MG_dep_mbs/k_mbs;
            MG_CONS_dep_kmeans_001=MG_CONS_dep_kmeans_001+current_MG_dep_kmeans_001/k_kmeans;
            MG_CONS_dep_mbs_001=MG_CONS_dep_mbs_001+current_MG_dep_mbs_001/k_mbs;

            MG_CONS_Time_kmeans=MG_CONS_Time_kmeans+current_MG_Time_kmeans/k_kmeans;
            MG_CONS_Time_mbs=MG_CONS_Time_mbs+current_MG_Time_mbs/k_mbs;
            MG_CONS_Time_kmeans_001=MG_CONS_Time_kmeans_001+current_MG_Time_kmeans_001/k_kmeans;
            MG_CONS_Time_mbs_001=MG_CONS_Time_mbs_001+current_MG_Time_mbs_001/k_mbs;

            MG_CONS_Time_kmeans_total=MG_CONS_Time_kmeans_total+current_MG_Time_kmeans;
            MG_CONS_Time_mbs_total=MG_CONS_Time_mbs_total+current_MG_Time_mbs;
            MG_CONS_Time_kmeans_001_total=MG_CONS_Time_kmeans_001_total+current_MG_Time_kmeans_001;
            MG_CONS_Time_mbs_001_total=MG_CONS_Time_mbs_001_total+current_MG_Time_mbs_001;

        end
        distribution_dep_kmeans=distribution_dep_kmeans+MG_CONS_dep_kmeans/groupLoop;
        distribution_dep_mbs=distribution_dep_mbs+MG_CONS_dep_mbs/groupLoop;
        distribution_dep_kmeans_001=distribution_dep_kmeans_001+MG_CONS_dep_kmeans_001/groupLoop;
        distribution_dep_mbs_001=distribution_dep_mbs_001+MG_CONS_dep_mbs_001/groupLoop;

        distribution_Time_kmeans=distribution_Time_kmeans+MG_CONS_Time_kmeans/groupLoop;
        distribution_Time_mbs=distribution_Time_mbs+MG_CONS_Time_mbs/groupLoop;
        distribution_Time_kmeans_001=distribution_Time_kmeans_001+MG_CONS_Time_kmeans_001/groupLoop;
        distribution_Time_mbs_001=distribution_Time_mbs_001+MG_CONS_Time_mbs_001/groupLoop;

        distribution_Time_kmeans_total=distribution_Time_kmeans_total+MG_CONS_Time_kmeans_total/groupLoop;
        distribution_Time_mbs_total=distribution_Time_mbs_total+MG_CONS_Time_mbs_total/groupLoop;
        distribution_Time_kmeans_001_total=distribution_Time_kmeans_001_total+MG_CONS_Time_kmeans_001_total/groupLoop;
        distribution_Time_mbs_001_total=distribution_Time_mbs_001_total+MG_CONS_Time_mbs_001_total/groupLoop;
    end 
     MG_DEP(index+1) = distribution_dep_kmeans/distributionLoop;
     UNIFORM_DEP(index+1) = distribution_dep_mbs/distributionLoop;
     MG_DEP_001(index+1) = distribution_dep_kmeans_001/distributionLoop;
     UNIFORM_DEP_001(index+1) = distribution_dep_mbs_001/distributionLoop;

     MG_TIME(index+1) = distribution_Time_kmeans/distributionLoop;
     UNIFORM_TIME(index+1) = distribution_Time_mbs/distributionLoop;
     MG_TIME_001(index+1) = distribution_Time_kmeans_001/distributionLoop;
     UNIFORM_TIME_001(index+1) = distribution_Time_mbs_001/distributionLoop;

     MG_TIME_TOTAL(index+1) = distribution_Time_kmeans_total/distributionLoop;
     UNIFORM_TIME_TOTAL(index+1) = distribution_Time_mbs_total/distributionLoop;
     MG_TIME_001_TOTAL(index+1) = distribution_Time_kmeans_001_total/distributionLoop;
     UNIFORM_TIME_001_TOTAL(index+1) = distribution_Time_mbs_001_total/distributionLoop;
end

delete(gcp);

figure(1); 
plot(loop_times,MG_DEP,'Color','k','LineStyle',':','Marker','o','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);
hold on;
plot(loop_times,UNIFORM_DEP,'Color','k','LineStyle','-','Marker','o','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);
plot(loop_times,MG_DEP_001,'Color','k','LineStyle',':','Marker','square','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);
plot(loop_times,UNIFORM_DEP_001,'Color','k','LineStyle','-','Marker','square','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);

ylim([0.3 0.44]);
yticks(0.3:0.02:0.44);
xlabel('Location estimation uncertainty of Willie, $\varepsilon_w^{2}$','Interpreter','latex');
ylabel('Average detection error probability');
legend('K-means++, $\epsilon=0.1$','SMBSP, $\epsilon=0.1$','K-means++, $\epsilon=0.2$','SMBSP, $\epsilon=0.2$','Location','best','Interpreter','latex');
box on;
grid on;

figure(2); 
plot(loop_times,MG_TIME,'Color','k','LineStyle',':','Marker','o','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);
hold on;
plot(loop_times,UNIFORM_TIME,'Color','k','LineStyle','-','Marker','o','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);
plot(loop_times,MG_TIME_001,'Color','k','LineStyle',':','Marker','square','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);
plot(loop_times,UNIFORM_TIME_001,'Color','k','LineStyle','-','Marker','square','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);

%ylim([0.3 0.5]);
%yticks(0.3:0.04:0.5);
xlabel('Location estimation uncertainty of Willie, $\varepsilon_w^{2}$','Interpreter','latex');
ylabel('Average transmission time (ms)');
legend('K-means++, $\epsilon=0.1$','SMBSP, $\epsilon=0.1$','K-means++, $\epsilon=0.2$','SMBSP, $\epsilon=0.2$','Location','best','Interpreter','latex');
box on;
grid on;

figure(3); 
plot(loop_times,MG_TIME_TOTAL,'Color','k','LineStyle',':','Marker','o','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);
hold on;
plot(loop_times,UNIFORM_TIME_TOTAL,'Color','k','LineStyle','-','Marker','o','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);
plot(loop_times,MG_TIME_001_TOTAL,'Color','k','LineStyle',':','Marker','square','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);
plot(loop_times,UNIFORM_TIME_001_TOTAL,'Color','k','LineStyle','-','Marker','square','MarkerFaceColor','none','MarkerSize',6,'MarkerIndices',1:1:length(loop_times),'linewidth',1);

%ylim([0.3 0.5]);
%yticks(0.3:0.04:0.5);
xlabel('Location estimation uncertainty of Willie, $\varepsilon_w^{2}$','Interpreter','latex');
ylabel('Total transmission time (ms)');
legend('K-means++, $\epsilon=0.1$','SMBSP, $\epsilon=0.1$','K-means++, $\epsilon=0.2$','SMBSP, $\epsilon=0.2$','Location','best','Interpreter','latex');
box on;
grid on;