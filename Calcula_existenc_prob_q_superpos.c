#include "mex.h"
#include "matrix.h"
#include "matrixutils.h"
#include "statutils.h"

#include <math.h>
//Subrutina para el problema de múltiple detección con etiquetas
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //Definicion de variables de entrada
    double  *Matriz_medidas, *Pos_sensores,  *SNR0, *d0, *Umbral_sensor,*Pos_med_blan;
    double *target_existence_prior_prob,*Existence_matrix,*indicador_blancos_act,*indicador_pos_med;
    //Definición de variables intermedias
    int j, k,  q, num_sensor;
    int ndim_1, Nsensores, Nsensores_cluster, Nblancos;
    double   *cdf;
    int *sensor_cluster, *num_sensor_cluster;

    mwSize dims_1[2];
    
    double d0_cuadrado, resta_x, resta_y, dist_cuadrado_sensor_q, Likelihood, SNR_sens, SNR_sens_sqrt;
    double u1, uj, peso_total,  peso_log_max;
    int ir, jr,i_a;
    int Nparticulas;
    
    double existence_prior_j;
    double hip_act;
    double *likelihood_existence;
    
    int Nexistence;
    //Definición de variables de salida
    double *existence_prob_q, *indicador_blancos_pred, *pesos_existence_pred,*hipotesis_particulas,*Nhipotesis,*Part_hipotesis,*hipotesis;
    
    //Inputs
    Matriz_medidas=mxGetPr(prhs[0]);
    Pos_sensores=mxGetPr(prhs[1]);
    SNR0=mxGetPr(prhs[2]);
    d0=mxGetPr(prhs[3]);
    Umbral_sensor=mxGetPr(prhs[4]);
    Pos_med_blan=mxGetPr(prhs[5]);
    target_existence_prior_prob=mxGetPr(prhs[6]);
    Existence_matrix=mxGetPr(prhs[7]);
    indicador_blancos_act=mxGetPr(prhs[8]);
    indicador_pos_med=mxGetPr(prhs[9]);
    
    
    //Dimensiones vectores
    Nsensores=(int) mxGetN(prhs[0]);
    Nblancos=(int) mxGetM(prhs[6]);
    Nexistence=(int) mxGetN(prhs[7]);
    Nparticulas=(int) mxGetN(prhs[8]);
    
    //Outputs
    ndim_1=2;
    dims_1[0] =1;
    dims_1[1] = Nexistence;
    plhs[0]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    existence_prob_q = mxGetPr(plhs[0]);
    
    ndim_1=2;
    dims_1[0] =Nblancos;
    dims_1[1] = Nparticulas;
    plhs[1]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    indicador_blancos_pred = mxGetPr(plhs[1]);
    
    ndim_1=2;
    dims_1[0] =1;
    dims_1[1] = Nparticulas;
    plhs[2]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    pesos_existence_pred = mxGetPr(plhs[2]);
    
    ndim_1=2;
    dims_1[0] =1;
    dims_1[1] = Nparticulas;
    plhs[3]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    hipotesis_particulas = mxGetPr(plhs[3]);
    
    ndim_1=2;
    dims_1[0] =1;
    dims_1[1] = 1;
    plhs[4]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    Nhipotesis = mxGetPr(plhs[4]);
    
    
    ndim_1=2;
    dims_1[0] =1;
    dims_1[1] = Nexistence+1;
    plhs[5]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    Part_hipotesis = mxGetPr(plhs[5]);
    
    ndim_1=2;
    dims_1[0] =1;
    dims_1[1] = Nexistence;
    plhs[6]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    hipotesis = mxGetPr(plhs[6]);
    
    
    //Reserva de memoria
    
    
    cdf=(double *)mxCalloc(Nexistence, sizeof(double));
    sensor_cluster=(int *)mxCalloc(Nsensores, sizeof(int));
    num_sensor_cluster=(int *)mxCalloc(Nsensores, sizeof(int));
    
    likelihood_existence=(double *)mxCalloc(Nexistence, sizeof(double));
    
    //Operaciones
    d0_cuadrado=d0[0]*d0[0];
    
    
    
    
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
    
    
    
    
    for (j=0;j<Nexistence;j++){
        existence_prior_j=1;
        
        for (k=0;k<Nblancos;k++){
            if(Existence_matrix[k+Nblancos*j]==1){
                existence_prior_j=existence_prior_j*target_existence_prior_prob[k];
            }else{
                existence_prior_j=existence_prior_j*(1-target_existence_prior_prob[k]);
                
            }
            
        }
        
        Likelihood=0;
        
        for (q=0;q<Nsensores_cluster;q++){
            num_sensor=num_sensor_cluster[q];
            //Calculamos la likelihood para este sensor
            SNR_sens=0;
            
            //Ahora se calcula la SNR de los blancos adicionales
            for (i_a=0;i_a<Nblancos;i_a++){
                if(Existence_matrix[i_a+Nblancos*j]==1){
                    resta_x=Pos_med_blan[2*i_a]-Pos_sensores[2*num_sensor];
                    resta_y=Pos_med_blan[2*i_a+1]-Pos_sensores[2*num_sensor+1];
                    dist_cuadrado_sensor_q=resta_x*resta_x+resta_y*resta_y;
                    if(dist_cuadrado_sensor_q>d0_cuadrado){
                        SNR_sens=SNR_sens+SNR0[0]*d0_cuadrado/dist_cuadrado_sensor_q;
                    }else{
                        SNR_sens=SNR_sens+SNR0[0];
                    }
                }
            }
            
            //Quito la sqrt
            //SNR_sens_sqrt=sqrt(SNR_sens);
            SNR_sens_sqrt=SNR_sens;
            Likelihood=Likelihood-(Matriz_medidas[num_sensor]-SNR_sens_sqrt)*(Matriz_medidas[num_sensor]-SNR_sens_sqrt)/2;
        }
        likelihood_existence[j]=Likelihood;
        
        existence_prob_q[j]=existence_prior_j;
        
    }
    
    //Se pasa a valor natural la likelihood
    peso_log_max=likelihood_existence[0];
    for (j=0;j<Nexistence;j++){
        if(likelihood_existence[j]>peso_log_max){
            peso_log_max=likelihood_existence[j];
        }
    }
    for (j=0;j<Nexistence;j++){
        likelihood_existence[j]=exp(-peso_log_max+likelihood_existence[j]);
    }
    //Se normaliza
    
    for (j=0;j<Nexistence;j++){
        existence_prob_q[j]=existence_prob_q[j]*likelihood_existence[j];
    }
    
    peso_total=0;
    for (j=0;j<Nexistence;j++){
        peso_total=peso_total+existence_prob_q[j];
    }
    
    for (j=0;j<Nexistence;j++){
        existence_prob_q[j]=existence_prob_q[j]/peso_total;
    }
    
    
    //Muestreo de la existence_prob_j (Necesito Nparticulas muestras)
    
    cdf[0]=existence_prob_q[0];
    //Este bucle empieza en 1
    for (j=1;j<Nexistence;j++){
        cdf[j]=cdf[j-1]+existence_prob_q[j];
    }
    
    u1=unirnd(0, 1)/Nparticulas;
    ir=0;
    
    for (jr=0;jr<Nparticulas;jr++){
        uj=u1+((double)jr)/Nparticulas;
        
        while(uj>cdf[ir]){
            ir=ir+1;
        }
        //Se hace el sampling
        for (k=0;k<Nblancos;k++){
            indicador_blancos_pred[k+Nblancos*jr]=Existence_matrix[k+Nblancos*ir];
        }
        pesos_existence_pred[jr]=likelihood_existence[ir];
        hipotesis_particulas[jr]=ir;
    }
    
    
    //Rellendo Nhipotesis
    Nhipotesis[0]=1;
    hip_act=hipotesis_particulas[0];
    hipotesis[0]=hip_act;
    Part_hipotesis[0]=0;

    for (jr=0;jr<Nparticulas;jr++){
        
        if(hipotesis_particulas[jr]>hip_act){
            Part_hipotesis[(int)Nhipotesis[0]]=jr;
            
            hip_act=hipotesis_particulas[jr];
            hipotesis[(int)Nhipotesis[0]]=hipotesis_particulas[jr];
            Nhipotesis[0]=Nhipotesis[0]+1;
            
            
        }
    }
    Part_hipotesis[(int)Nhipotesis[0]]=Nparticulas;
    
 
}

