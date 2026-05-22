#include "mex.h"
#include "matrix.h"
#include "matrixutils.h"
#include "statutils.h"

#include <math.h>
//Subrutina para el problema de múltiple detección con etiquetas
//Basada en subrutina_lmb_pos_pp_sin_rao2
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //Definicion de variables de entrada
    double *particulas_act, *Matriz_medidas, *Pos_sensores,  *phi, *epsilon, *Umbral_sensor, *t_vida, *Pos_med_blan, *T;
    double *chol_Q,*A,*indicador_blancos_act,*indicador_blancos_pred,*indicador_pos_med;
    double *hipotesis_particulas,*Nhipotesis,*Part_hipotesis,*hipotesis,*Existence_matrix,*vector_aleatorio_global;
    double *indicador_b,*means_b,*chol_P_ini;
    //Definición de variables intermedias
    int j, k,  q, num_sensor;
    int ndim_1, Nsensores, Nsensores_cluster, Nblancos;

    mwSize dims_1[2]; //This was the key line to change to make old mex files work!!

    double  *pesos_pred_inter, *cdf, *particulas_inter;
    int *sensor_cluster, *num_sensor_cluster;
    double  ganancia_vel_act[4], chol_act[4], vector_aleatorio[4], posicionp_act[2], posicionp_pred[2], posicionv_act[2], innov_pos[2],vector_aleatorio_mult[4];
    double termino4[2], vel_fin[2],particula_k_act[4],particula_k_pred[4];
    double resta_x, resta_y, dist_cuadrado_sensor_q, Likelihood, SNR_sens, SNR_sens_sqrt;
    double u1, uj, peso_total,  peso_log_max;
    int ir, jr,i_a;
    int Nparticulas,  t_vida_i,Nparticulas_h;
    int h;
    int Nexistence,hipotesis_h,Nparticulas_k_1;
    double suma_bias;
    int Nhipotesis_int,index_aux_b;
    
    //Definición de variables de salida
    double *particulas_pred, *Bias,*index_parent;;
    
    //Inputs
    particulas_act=mxGetPr(prhs[0]);
    Matriz_medidas=mxGetPr(prhs[1]);
    Pos_sensores=mxGetPr(prhs[2]);
    phi=mxGetPr(prhs[3]);
    epsilon=mxGetPr(prhs[4]);
    Umbral_sensor=mxGetPr(prhs[5]);
    t_vida=mxGetPr(prhs[6]);
    Pos_med_blan=mxGetPr(prhs[7]);
    T=mxGetPr(prhs[8]);
    chol_Q=mxGetPr(prhs[9]);
    A=mxGetPr(prhs[10]); // A es la matriz de transicion (modelo dinamico)
    indicador_blancos_act=mxGetPr(prhs[11]);
    indicador_blancos_pred=mxGetPr(prhs[12]);
    indicador_pos_med=mxGetPr(prhs[13]);
    hipotesis_particulas=mxGetPr(prhs[14]);
    Nhipotesis=mxGetPr(prhs[15]);
    Part_hipotesis=mxGetPr(prhs[16]);
    hipotesis=mxGetPr(prhs[17]);
    Existence_matrix=mxGetPr(prhs[18]);
    vector_aleatorio_global=mxGetPr(prhs[19]);
    indicador_b=mxGetPr(prhs[20]);
    means_b=mxGetPr(prhs[21]);
    chol_P_ini=mxGetPr(prhs[22]);
    
    
    
    //Dimensiones vectores
    Nparticulas=mxGetN(prhs[0]);
    Nsensores=mxGetN(prhs[1]);
    Nblancos=mxGetM(prhs[0])/4;
    Nexistence=mxGetN(prhs[18]);
    
    t_vida_i=(int)(t_vida[0]-1);//Se destaca que se resta 1
    //Outputs
    ndim_1=2;
    dims_1[0] =4*Nblancos;
    dims_1[1] = Nparticulas;
    plhs[0]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    particulas_pred = mxGetPr(plhs[0]);
    
    ndim_1=2;
    dims_1[0]=Nparticulas;
    dims_1[1]=Nblancos;
    plhs[1]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    Bias=mxGetPr(plhs[1]);
    
    ndim_1=2;
    dims_1[0]=Nparticulas;
    dims_1[1]=Nblancos;
    plhs[2]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    index_parent=mxGetPr(plhs[2]);
    
    
    //Reserva de memoria
    
    pesos_pred_inter= (double *) mxCalloc(Nparticulas, sizeof(double));
    cdf=(double *)mxCalloc(Nparticulas, sizeof(double));
    sensor_cluster=(int *)mxCalloc(Nsensores, sizeof(int));
    num_sensor_cluster=(int *)mxCalloc(Nsensores, sizeof(int));
    particulas_inter=(double *)mxCalloc(4*Nparticulas, sizeof(double));
    
    //Operaciones
    
    
    //Se calculan los sensores que pertenecen al cluster
    for (k=0;k<Nblancos;k++){
        if(indicador_pos_med[k]==1){
            for(q=0;q<Nsensores;q++){
                resta_x=Pos_med_blan[2*k]-Pos_sensores[2*q];
                resta_y=Pos_med_blan[2*k+1]-Pos_sensores[2*q+1];
                dist_cuadrado_sensor_q=resta_x*resta_x+resta_y*resta_y;
                if(dist_cuadrado_sensor_q<Umbral_sensor[0]){
                    sensor_cluster[q]=1;
                }
            }
        }
    }
    
    
    Nsensores_cluster=0;
    //Se recorre sensor_cluster para ver el numero de sensores que tenemos en el cluster y sus id
    for (q=0;q<Nsensores;q++){
        if(sensor_cluster[q]==1){
            num_sensor_cluster[Nsensores_cluster]=q;
            Nsensores_cluster=Nsensores_cluster+1;
        }
    }
    
    
    Nhipotesis_int=(int)Nhipotesis[0];
    
    for(h=0;h<Nhipotesis[0];h++){
        
        hipotesis_h=(int)hipotesis[h];
        
        for (k=0;k<Nblancos;k++){
            Nparticulas_k_1=0;
            if(Existence_matrix[k+Nblancos*hipotesis_h]){                   
                    for (j=0;j<Nparticulas;j++){
                        if(indicador_blancos_act[k+Nblancos*j]){
                            
                            Nparticulas_k_1=Nparticulas_k_1+1;
                            
                            // We take the random vector created in Matlab and passed as an input
                            vector_aleatorio[0]=vector_aleatorio_global[h+Nhipotesis_int*k+Nhipotesis_int*Nblancos*j+Nhipotesis_int*Nblancos*Nparticulas*0];
                            vector_aleatorio[1]=vector_aleatorio_global[h+Nhipotesis_int*k+Nhipotesis_int*Nblancos*j+Nhipotesis_int*Nblancos*Nparticulas*1];
                            vector_aleatorio[2]=vector_aleatorio_global[h+Nhipotesis_int*k+Nhipotesis_int*Nblancos*j+Nhipotesis_int*Nblancos*Nparticulas*2];
                            vector_aleatorio[3]=vector_aleatorio_global[h+Nhipotesis_int*k+Nhipotesis_int*Nblancos*j+Nhipotesis_int*Nblancos*Nparticulas*3];
                         
                            if(indicador_b[k]==1){
                                //Se trata de un blanco nuevo
                                //Hay que ver que blanco es y ponerlo con particula_k_pred
                                index_aux_b=0;
                                for(q=0;q<k;q++){
                                    if(indicador_b[q]==1){
                                        index_aux_b=index_aux_b+1;
                                    }
                                }
                                
                                particula_k_pred[0]=means_b[0+4*index_aux_b];
                                particula_k_pred[1]=means_b[1+4*index_aux_b];
                                particula_k_pred[2]=means_b[2+4*index_aux_b];
                                particula_k_pred[3]=means_b[3+4*index_aux_b];
                                        
                                matmult(chol_P_ini,4,4,vector_aleatorio,4, 1,vector_aleatorio_mult);
        
                                
                            }else{
                                //Se trata de un blanco antiguo
                                
                                posicionp_act[0]=particulas_act[4*k+4*Nblancos*j];
                                posicionp_act[1]=particulas_act[4*k+2+4*Nblancos*j];
                                posicionv_act[0]=particulas_act[4*k+1+4*Nblancos*j];
                                posicionv_act[1]=particulas_act[4*k+3+4*Nblancos*j];
                                
                                //No tiene Rao Blackwellisation
                                
                                particula_k_act[0]=posicionp_act[0];
                                particula_k_act[1]=posicionv_act[0];
                                particula_k_act[2]=posicionp_act[1];
                                particula_k_act[3]=posicionv_act[1];
                                
                                
                                
                                
                                matmult(chol_Q,4,4,vector_aleatorio,4, 1,vector_aleatorio_mult);
                                matmult(A,4,4,particula_k_act,4, 1,particula_k_pred);
                                
                            }
                            
                            particula_k_pred[0]=particula_k_pred[0]+vector_aleatorio_mult[0];
                            particula_k_pred[1]=particula_k_pred[1]+vector_aleatorio_mult[1];
                            particula_k_pred[2]=particula_k_pred[2]+vector_aleatorio_mult[2];
                            particula_k_pred[3]=particula_k_pred[3]+vector_aleatorio_mult[3];
                            
                            
                            posicionp_pred[0]=particula_k_pred[0];
                            posicionp_pred[1]=particula_k_pred[2];
                            vel_fin[0]=particula_k_pred[1];
                            vel_fin[1]=particula_k_pred[3];
                            
                            
                            //Se calcula la likelihood en escala logaritmica
                            Likelihood=0;
                            
                            for (q=0;q<Nsensores_cluster;q++){
                                num_sensor=num_sensor_cluster[q];
                                //Calculamos la likelihood para este sensor
                                SNR_sens=0;
                                //Primero para el blanco actual
                                
                                resta_x=posicionp_pred[0]-Pos_sensores[2*num_sensor];
                                resta_y=posicionp_pred[1]-Pos_sensores[2*num_sensor+1];
                                dist_cuadrado_sensor_q=resta_x*resta_x+resta_y*resta_y;
                               
                                
                                SNR_sens=SNR_sens+phi[0]/(dist_cuadrado_sensor_q+epsilon[0]);

                                
//                                 if(dist_cuadrado_sensor_q>d0_cuadrado){
//                                     SNR_sens=SNR_sens+SNR0[0]*d0_cuadrado/dist_cuadrado_sensor_q;
//                                 }else{
//                                     SNR_sens=SNR_sens+SNR0[0];
//                                 }
                                //Ahora se calcula la SNR de los blancos adicionales
                                for (i_a=0;i_a<Nblancos;i_a++){
                                    
                                    if(Existence_matrix[i_a+Nblancos*hipotesis_h]){//Tiene que pertenecer a la hipotesis (todas estas particulas tienen el mismo indicador_blancos_pred
                                        if(i_a!=k){
                                            
                                            resta_x=Pos_med_blan[2*i_a]-Pos_sensores[2*num_sensor];
                                            resta_y=Pos_med_blan[2*i_a+1]-Pos_sensores[2*num_sensor+1];
                                            dist_cuadrado_sensor_q=resta_x*resta_x+resta_y*resta_y;


                                            SNR_sens=SNR_sens+phi[0]/(dist_cuadrado_sensor_q+epsilon[0]);


//                                             if(dist_cuadrado_sensor_q>d0_cuadrado){
//                                                 SNR_sens=SNR_sens+SNR0[0]*d0_cuadrado/dist_cuadrado_sensor_q;
//                                             }else{
//                                                 SNR_sens=SNR_sens+SNR0[0];
//                                             }
                                        }
                                    }
                                }
                                
                                //Quito sqrt
                                //SNR_sens_sqrt=sqrt(SNR_sens);
                                SNR_sens_sqrt=SNR_sens;
                                
                                Likelihood=Likelihood-(Matriz_medidas[num_sensor]-SNR_sens_sqrt)*(Matriz_medidas[num_sensor]-SNR_sens_sqrt)/2;
                            }
                            pesos_pred_inter[j]=Likelihood;
                            particulas_inter[4*j]=posicionp_pred[0];
                            particulas_inter[4*j+1]=vel_fin[0];
                            particulas_inter[4*j+2]=posicionp_pred[1];
                            particulas_inter[4*j+3]=vel_fin[1];
                        }else{
                            
                            pesos_pred_inter[j]=0; //Hay que tener esto en cuenta en la normalizacion de los pesos
                        }
                    }
                
                
                //Hay que pasar a valor natural los pesos de las particulas
                
                //Se coge la primera particula con peso menor que cero
                peso_log_max=0;
                j=0;
                while(peso_log_max==0){
                    if(pesos_pred_inter[j]<0){
                        peso_log_max=pesos_pred_inter[j];
                    }
                    j=j+1;
                }
                
                for (j=0;j<Nparticulas;j++){
                    if(pesos_pred_inter[j]<0){
                        if(pesos_pred_inter[j]>peso_log_max){
                            peso_log_max=pesos_pred_inter[j];
                        }
                    }
                }
                for (j=0;j<Nparticulas;j++){
                    if(pesos_pred_inter[j]<0){
                        pesos_pred_inter[j]=exp(-peso_log_max+pesos_pred_inter[j]);
                    }
                }
                
                
                //Hay que normalizar los pesos de las particulas
                
                peso_total=0;
                for (j=0;j<Nparticulas;j++){
                    peso_total=peso_total+pesos_pred_inter[j];
                }
                for(j=0;j<Nparticulas;j++){
                    pesos_pred_inter[j]=pesos_pred_inter[j]/peso_total;
                }
                
                
                //Se muestrean los índices del filtro auxiliar
                //Se calcula la cdf
                
                cdf[0]=pesos_pred_inter[0];
                //Este bucle empieza en 1
                for (j=1;j<Nparticulas;j++){
                    cdf[j]=cdf[j-1]+pesos_pred_inter[j];
                }
                
                //Se inicia con la variable aleatoria correspondiente
                
                Nparticulas_h=Part_hipotesis[h+1]-Part_hipotesis[h]; //Este va indexado con h (no con hipotesis_h)
                
                u1=unirnd(0, 1)/Nparticulas_h;
                ir=0;
                
                suma_bias=0;
                
                for (jr=Part_hipotesis[h];jr<Part_hipotesis[h+1];jr++){
                    
                    uj=u1+((double)(jr-Part_hipotesis[h]))/Nparticulas_h;
                    
                    while(uj>cdf[ir]){
                        ir=ir+1;
                    }
                    //Ahora hay que sacar las nuevas particulas
                    particulas_pred[4*k+4*Nblancos*jr]=particulas_inter[4*ir];
                    particulas_pred[4*k+1+4*Nblancos*jr]=particulas_inter[4*ir+1];
                    particulas_pred[4*k+2+4*Nblancos*jr]=particulas_inter[4*ir+2];
                    particulas_pred[4*k+3+4*Nblancos*jr]=particulas_inter[4*ir+3];
                    Bias[jr+Nparticulas*k]=pesos_pred_inter[ir];//En pp hay que multiplicar pues es un Bias por hipótesis no por blanco
                    index_parent[jr+Nparticulas*k]=ir;
                    
                    suma_bias=suma_bias+Bias[jr+Nparticulas*k];
                    
                }
                
                for (jr=Part_hipotesis[h];jr<Part_hipotesis[h+1];jr++){
                    //Bias[jr+Nparticulas*k]=Bias[jr+Nparticulas*k]/suma_bias;
                    
                    Bias[jr+Nparticulas*k]=Bias[jr+Nparticulas*k]*Nparticulas_k_1;
                }
                
                
            }
        }
    }
}
