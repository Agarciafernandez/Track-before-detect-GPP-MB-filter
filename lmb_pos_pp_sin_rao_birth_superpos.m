function [weights_pred,particulas_pred,indicador_blancos_pred,index_parent,mean_particles_pred,Pos_med_blan,indicador_pos_med]=lmb_pos_pp_sin_rao_birth_superpos(particles_act,...
    Matrix_measurements,Pos_sensors,phi,epsilon,Threshold_sensor,k,T,weights_act,chol_Q,A,indicator_targets_act,Pos_med_blan,...
    p_s,Existence_matrix,indicador_pos_med,s_rand2,Ncom_b,means_b,chol_P_ini,p_b,Pos_med_blan_b)


%This filter runs the sampling and particle weight calculation in 
%[1] A. F. Garcia-Fernandez, "Track-before-detect labeled multi-Bernoulli particle filter with label switching," 
% in IEEE Transactions on Aerospace and Electronic Systems, vol. 52, no. 5, pp. 2123-2138, October 2016

%See Algorithm 2 in [1].

%Maximum number of targets
Ntargets_max=size(particles_act,1)/4;

%We calculate the predicted states of the targets
mean_particles_pred=kron(eye(Ntargets_max),A)*particles_act;

%We calculate the predicted probabilities of existence
Nparticles=size(indicator_targets_act,2);
target_existence_prior_prob=p_s*sum(indicator_targets_act,2)/Nparticles;
%We add target births (looking for Ncom_b free spaces, first in the existence variable and then in Pos_med_blan)


index_b=find(target_existence_prior_prob==0);
index_b=index_b(1:Ncom_b)';
indicador_b=zeros(size(indicador_pos_med));
indicador_b(index_b)=1;
target_existence_prior_prob(index_b)=p_b;
aux_index_b=[2*index_b-1;2*index_b];
aux_index_b=aux_index_b(:);

indicador_pos_med_old=indicador_pos_med;
Pos_med_blan(aux_index_b)=Pos_med_blan_b;
indicador_pos_med(index_b)=1;



%We sample set of labels using the predicted positions for the targets
%using (34), and the different label sets with its number of particles.
%Nhypotheses is the different number of sets of labels
[~,indicador_blancos_pred,weights_existence_pred,hypotheses_particles,Nhypotheses,Part_hypothesis,hypothesis]=Calcula_existenc_prob_q_superpos_update(Matrix_measurements,Pos_sensors,phi,epsilon,Threshold_sensor,...
    Pos_med_blan,target_existence_prior_prob,Existence_matrix,indicator_targets_act,indicador_pos_med);




%Solo considero los sensores proximos a las zonas donde muestreo (No hace
%falta tocar Pos_med_blan)
indicador_pos_med=double(sum(indicador_blancos_pred,2)>1);
%Modifico indicador_blancos_act para meter los target birth
indicador_pos_med_birth=and(indicador_pos_med==1,indicador_pos_med_old==0);

indicator_targets_act(indicador_pos_med_birth,:)=1;


%Auxiliary routine. We sample the particles and calculate the bias term,
%also performing resampling per target, see Algorithm 2. We draw the random
%vector outside C to have better understanding.
random_vector_global=randn(s_rand2,1,Nhypotheses*Ntargets_max*length(weights_act)*4);
[particulas_pred,Bias,index_parent]=subrutina_lmb_pos_pp_sin_rao_birth_superpos_update(particles_act,Matrix_measurements,Pos_sensors,phi,epsilon,Threshold_sensor,...
    k,Pos_med_blan,T,chol_Q,A,indicator_targets_act,indicador_blancos_pred,indicador_pos_med,hypotheses_particles,Nhypotheses,...
    Part_hypothesis,hypothesis,Existence_matrix,random_vector_global,indicador_b,means_b,chol_P_ini);


%We compute the particle weight, see (38) in 

weights_pred=Peso_lmb_pos_pp_superpos_update(particulas_pred,Matrix_measurements,Pos_sensors,...
    phi,epsilon,Threshold_sensor,Pos_med_blan,T,Bias,indicador_blancos_pred,weights_existence_pred,indicador_pos_med);


