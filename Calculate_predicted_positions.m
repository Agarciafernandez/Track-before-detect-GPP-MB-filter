function [Pos_med_blan,indicador_pos_med]=Calculate_predicted_positions(pesos_act,particulas_act,indicador_blancos_act,T)
%This function groups the targets in clusters based on the predicted
%positions of the targets that have some alive particles at the previous
%time steps

Pos_med_blan=zeros(2,size(indicador_blancos_act,1));

indicador_pos_med=double(sum(indicador_blancos_act,2)>0);

A=[1 T 0 0; 0 0 1 T];

blancos_activos_ini=find(sum(indicador_blancos_act,2)>0);
Nblancos=length(blancos_activos_ini);



if(Nblancos>0)
    
    for i=1:Nblancos
        num_blanco=blancos_activos_ini(i);
        ind_part_blanco=find(indicador_blancos_act(num_blanco,:));
        pesos_con_blanco=pesos_act(ind_part_blanco);
        %Hay que normalizar los pesos de cada blanco
        pesos_con_blanco=pesos_con_blanco/sum(pesos_con_blanco);
        part_con_blanco=particulas_act(4*num_blanco-3:4*num_blanco, ind_part_blanco);
        
        Pos_vel_t_1=part_con_blanco*pesos_con_blanco';
        Pos_med_blan(:,num_blanco)=A*Pos_vel_t_1;
        
    end
    
    
end