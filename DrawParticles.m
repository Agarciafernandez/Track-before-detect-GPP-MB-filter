function DrawParticles(particles_act,indicator_targets_act,X_multi_k,Area,k,t_birth,t_death)

%This function plots the particles of each target and the ground truth

clf
axis([0 Area(1),0 Area(2)])
 
%  axis([X_multi_k(1)-5 X_multi_k(1)+5, X_multi_k(3)-5 X_multi_k(3)+5])

grid on
hold on

targets_active_ini=find(sum(indicator_targets_act,2)>0);
Ntargets=length(targets_active_ini);

for i=1:Ntargets
    num_target=targets_active_ini(i);
    ind_part_target=find(indicator_targets_act(num_target,:));
    plot(particles_act((num_target-1)*4+1,ind_part_target),particles_act((num_target-1)*4+3,ind_part_target),'o')
 
end


for i=1:size(X_multi_k,1)/4
    if(and(k>=t_birth(i),k<t_death(i)))

        plot(X_multi_k((i-1)*4+1),X_multi_k((i-1)*4+3),'r')
        text(X_multi_k((i-1)*4+1),X_multi_k((i-1)*4+3),num2str(i)) 
    end
end



