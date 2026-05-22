%Implementation of the generalised parallel partition multi-Bernoulli particle filter (without label switching) described in  
%[1] A. F. Garcia-Fernandez, "Track-before-detect labeled multi-Bernoulli particle filter with label switching," 
% in IEEE Transactions on Aerospace and Electronic Systems, vol. 52, no. 5, pp. 2123-2138, October 2016

%The Measurement model in this script is slightly different from [1] and is
%the one with superpositional sensors in [2]

%[2] E. S. Davies and Á. F. García-Fernández, "Information Exchange Track-Before-Detect Multi-Bernoulli Filter for Superpositional Sensors," in IEEE Transactions on Signal Processing, vol. 72, pp. 607-621, 2024,doi: 10.1109/TSP.2024.3349769

%This scenario considers the informative birth scenario. This version does
%not perform clustering.



%Author: A. F. Garcia-Fernandez

clear

addpath('IEMB files')
addpath('GOSPA')

rng(2,'twister')
%number of targets
numtruth = 5;

%create Sensor Network
Area=[120 120];
Cells=[12 12];
cell_width = Area(1)/Cells(1);

%initalise measurement vector
[sensor_pos , sensor_index]  = createSensorNetwork(Cells,Area);
sensor_num = length(sensor_pos);

%set state space dimension
n=4;

%Nearly constant velocity model
T=1;
F=[1 T 0 0;0 1 0 0; 0 0 1 T;0 0 0 1];
sigmaU=0.5;
Q=(sigmaU)^2*[T^3/3 T^2/2 0 0; T^2/2 T 0 0;0 0 T^3/3 T^2/2; 0 0 T^2/2 T];

%Set the measurement covariance
R=1*eye(Cells(1)*Cells(2));
p_s_constant=0.99;

%Number of steps
Nsteps=81;

nmc=100; % Number of Monte Carlo Loops

%We generate ground truth according to the model
[X_truth] = Trajectory_cross(Nsteps, F, numtruth, Q, Area);

% Makes sure the target trajectories are contained within the surveillence region
[X_truth,m_b,t_birth,t_death] = Trajectory_region(X_truth,numtruth,Nsteps);

%Round numbers for the birth model

m_b=[75,0,100,0;...
    20,0,105,0;...
    110,0,50,0;...
    110,0,85,0]';

P_b=diag([10 10 10 10]);
P_ini=repmat(P_b,1,1,numtruth);





%Define variables particle filter
Ntargets_max=15;
Nparticles=10000;
Pos_sensors=sensor_pos; %It has been checked that Pos_sensores=Calcula_Pos_sensores(Area,[cell_width,cell_width], 0,0,144); coincides with this line



%Measurement parameters used in [2] (instead of SNR0,d0)
phi_measurement = 500;
epsilon_measurement= 25;


Threshold_sensor=40^2; %The sensors whose distance to the mean positions in each cluster of targets is smaller than this threshold are considered in the update
Existence_matrix=CalculateExistenceMatrix(Ntargets_max); 
chol_Q=chol(Q)';
p_s=p_s_constant;
chol_P_ini=chol(P_b)';
means_b=m_b;
Ncom_b=size(m_b,2);
p_b= 0.0001;

Pos_med_blan_b=[means_b(1,:);means_b(3,:)];
Pos_med_blan_b=Pos_med_blan_b(:);

index_mult_existence=[1,cumprod(2*ones(1,Ntargets_max-1))];



Ntargets_max_est=8;
Perm_matrix=perms(0:Ntargets_max_est-1);
Perm_matrix=fliplr(Perm_matrix);

[s_rand1,s_rand2] = RandStream.create('mrg32k3a','NumStreams',3,'Seed',0);


%Variables for GOSPA
p=2;
alpha=2;


%Define GOSPA variables
mb_pf_squared_gospa_t_tot=zeros(1,Nsteps); %Total error
mb_pf_squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
mb_pf_squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
mb_pf_squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error

nmc=1

for q=1 :nmc

    tic

    z=measure(X_truth,t_birth,t_death,R,Nsteps,n,sensor_num,sensor_pos,sensor_index,Cells);

    %Particles and weights
    particles_act=zeros(4*Ntargets_max,Nparticles);
    weights_act=ones(1,Nparticles)/Nparticles;

    %Indicator to indicate if the particles contain the targets

    indicator_targets_act=zeros(Ntargets_max,Nparticles);

    for k=1:Nsteps

        heat=zeros(12,12);
        for i=1:144
            heat(sensor_index(1,i),sensor_index(2,i))= z(i,k);
        end
        Matrix_measurements= reshape(heat(:,:),(size(heat,2)*size(heat,1)),1)';


        [Pos_med_blan,indicador_pos_med]=Calculate_predicted_positions(weights_act,particles_act,indicator_targets_act,T); 
        Pos_med_blan=Pos_med_blan(:);     


        %We run the GPP method

        [weights_pred,particles_pred,indicator_targets_pred,index_parent_pred,mean_particulas_pred,Pos_med_blan,indicador_pos_med]=lmb_pos_pp_sin_rao_birth_superpos(particles_act,...
            Matrix_measurements,Pos_sensors,phi_measurement,epsilon_measurement,Threshold_sensor,k,T,weights_act,chol_Q,F,indicator_targets_act,Pos_med_blan,...
            p_s,Existence_matrix,indicador_pos_med,s_rand2,Ncom_b,means_b,chol_P_ini,p_b,Pos_med_blan_b);

        %Resample step
        [particles_pred_res,index_parent_pred_res,indicator_targets_pred_res]=resamplePF_lmb_index(particles_pred,weights_pred,index_parent_pred,indicator_targets_pred);
        
        
        weights_act=ones(1,Nparticles)/Nparticles;
        particles_act=particles_pred_res;
        indicator_targets_act=indicator_targets_pred_res;

        %Uncomment to plot particles
        % figure(1)
        % DrawParticles(particles_act,indicator_targets_act,X_truth(:,k),Area,k,t_birth,t_death)

        %We calculate the cardinality probability mass function
        [~,prob_cardinality]=Calcula_prob_hypothesis_lmb(particles_act,indicator_targets_act,Existence_matrix,index_mult_existence);

        %We estimate the cardinality and then the position of the targets
        %using kmeans
        [max_prob,cardinality_est]=max(prob_cardinality);
        cardinality_est=cardinality_est-1;
        indices_card=double(sum(indicator_targets_act,1)==cardinality_est);
        Nparticulas_card=sum(indices_card);

        particulas_act_card=Calcula_particulas_cardinality(particles_act,indicator_targets_act,cardinality_est,indices_card,Nparticulas_card);

        Perm_matrix_card=Perm_matrix(1:factorial(cardinality_est),1:cardinality_est);
        [particulas_kmeans,targets_estimate,indicador_cambio]=kmeans_mex3(particulas_act_card, ones(1,Nparticulas_card)/Nparticulas_card,Perm_matrix_card,50);
        X_estimate_pos= [targets_estimate(1:4:end)';targets_estimate(3:4:end)'];

        


        %We obtain the alive targets and compute GOSPA metric        
        alive_targets=and(t_birth<=k,t_death>=k); %get the index positions for detected targets
        X_truth_pos = reshape([alive_targets;and(false,alive_targets);alive_targets;and(false,alive_targets)],[],1);
        X_truth_pos = reshape(X_truth(X_truth_pos,k),2,[]);
        
        [d_gospa, ~, decomp_cost] = GOSPA(X_truth_pos, X_estimate_pos, p, cell_width/2, alpha);

        squared_gospa=d_gospa^2;
        gospa_loc=decomp_cost.localisation;
        gospa_mis=decomp_cost.missed;
        gospa_fal=decomp_cost.false;

        mb_pf_squared_gospa_t_tot(k)=mb_pf_squared_gospa_t_tot(k)+squared_gospa;
        mb_pf_squared_gospa_loc_t_tot(k)=mb_pf_squared_gospa_loc_t_tot(k)+gospa_loc;
        mb_pf_squared_gospa_false_t_tot(k)=mb_pf_squared_gospa_false_t_tot(k)+gospa_fal;
        mb_pf_squared_gospa_mis_t_tot(k)=mb_pf_squared_gospa_mis_t_tot(k)+gospa_mis;




    end
    t=toc;
    display(['Completed iteration number ', num2str(q),' time ', num2str(t), ' sec'])

end


%Root mean square GOSPA errors at each time step for IMB-UKF
mb_pf_rms_gospa_t=sqrt(mb_pf_squared_gospa_t_tot/nmc);
mb_pf_rms_gospa_loc_t=sqrt(mb_pf_squared_gospa_loc_t_tot/nmc);
mb_pf_rms_gospa_false_t=sqrt(mb_pf_squared_gospa_false_t_tot/nmc);
mb_pf_rms_gospa_mis_t=sqrt(mb_pf_squared_gospa_mis_t_tot/nmc);


%save('GPP_MB_PF_informative_results')


figure(1)
clf
plot(1:Nsteps,mb_pf_rms_gospa_t,'blue','Linewidth',1.3)
lim =axis;
ax=[0 Nsteps lim(3:4)]; %Ensures all figures have the same axes
axis(ax);
grid on
xlabel('Time step')
ylabel('RMS GOSPA error')
legend('Location','northwest')
set(gca,'FontSize',15)


figure(2)
clf
plot(1:Nsteps,mb_pf_rms_gospa_loc_t,'blue','Linewidth',1.3)
hold on
axis(ax);
grid on
xlabel('Time step')
ylabel('RMS GOSPA localisation error')
legend('Location','northwest')
set(gca,'FontSize',15)

figure(3)
clf
plot(1:Nsteps,mb_pf_rms_gospa_false_t,'blue','Linewidth',1.3)
axis(ax);
grid on
xlabel('Time step')
ylabel('RMS GOSPA false target error')
legend('Location','northwest')
set(gca,'FontSize',15)

figure(4)
clf
plot(1:Nsteps,mb_pf_rms_gospa_mis_t,'blue','Linewidth',1.3)
axis(ax);
grid on
xlabel('Time step')
ylabel('RMS GOSPA missed target error')
legend('Location','northwest')
set(gca,'FontSize',15)
