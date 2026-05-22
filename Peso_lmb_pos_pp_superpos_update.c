#include "mex.h"
#include "matrix.h"
#include "matrixutils.h"
#include <math.h>
//Subrutina para el problema de múltiple detección con etiquetas
//utiliza una likelihood para cada partícula
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //Definicion de variables de entrada
    double *particulas_pred, *Matriz_medidas, *Pos_sensores, *phi, *epsilon;
    double *Umbral_sensor,  *Pos_med_blan, *T, *Bias, *indicador_blancos_pred,*pesos_existence_pred,*indicador_pos_med;
    
    //Definicion de variables intermedias
    int  Nparticulas, ndim_1, i, j, k, q, Nsensores_cluster, Nblancos, num_sensor;
    int Nsensores;
    mwSize dims_1[2]; //This was the key line to change to make old mex files work!!
    double posicionp_pred[2];
    double dist_cuadrado_sensor_q, resta_x, resta_y, SNR_sens, SNR_sens_sqrt,Likelihood;
    int *sensor_cluster, *num_sensor_cluster;
    double peso_total,  peso_log_max,producto_bias_log;
    
    //Definición de variables de salida
    double *pesos_pred;
    
    //Inputs
    
    particulas_pred=mxGetPr(prhs[0]);
    Matriz_medidas=mxGetPr(prhs[1]);
    Pos_sensores=mxGetPr(prhs[2]);
    phi=mxGetPr(prhs[3]);
    epsilon=mxGetPr(prhs[4]);
    Umbral_sensor=mxGetPr(prhs[5]);
    Pos_med_blan=mxGetPr(prhs[6]);
    T=mxGetPr(prhs[7]);
    Bias=mxGetPr(prhs[8]);
    indicador_blancos_pred=mxGetPr(prhs[9]);
    pesos_existence_pred=mxGetPr(prhs[10]);
    indicador_pos_med=mxGetPr(prhs[11]);
    
    //Dimensiones vectores
    Nparticulas=mxGetN(prhs[0]);
    Nblancos=mxGetM(prhs[0])/4;
    Nsensores=mxGetN(prhs[1]);
    
    
    //Outputs
    ndim_1=2;
    dims_1[0] =1;
    dims_1[1] = Nparticulas;
    plhs[0]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    pesos_pred = mxGetPr(plhs[0]);
    
    
    //Reserva de memoria
    sensor_cluster=(int *)mxCalloc(Nsensores, sizeof(int));
    num_sensor_cluster=(int *)mxCalloc(Nsensores, sizeof(int));
    
    //Operaciones
    
    
    
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
    
    
    
    for (j=0;j<Nparticulas;j++){
        Likelihood=0;
        
        for (q=0;q<Nsensores_cluster;q++){
            num_sensor=num_sensor_cluster[q];
            //Calculamos la likelihood para este sensor
            SNR_sens=0;
            for (k=0;k<Nblancos;k++){
                if(indicador_blancos_pred[k+Nblancos*j]){
                    posicionp_pred[0]=particulas_pred[4*k+4*Nblancos*j];
                    posicionp_pred[1]=particulas_pred[4*k+2+4*Nblancos*j];
                    
                    resta_x=posicionp_pred[0]-Pos_sensores[2*num_sensor];
                    resta_y=posicionp_pred[1]-Pos_sensores[2*num_sensor+1];
                    dist_cuadrado_sensor_q=resta_x*resta_x+resta_y*resta_y;

                    SNR_sens=SNR_sens+phi[0]/(dist_cuadrado_sensor_q+epsilon[0]);


//                     if(dist_cuadrado_sensor_q>d0_cuadrado){
//                         SNR_sens=SNR_sens+SNR0[0]*d0_cuadrado/dist_cuadrado_sensor_q;
//                     }else{
//                         SNR_sens=SNR_sens+SNR0[0];
//                     }
                }
            }
            //Quito sqrt
            //SNR_sens_sqrt=sqrt(SNR_sens);
            SNR_sens_sqrt=SNR_sens;
            Likelihood=Likelihood-(Matriz_medidas[num_sensor]-SNR_sens_sqrt)*(Matriz_medidas[num_sensor]-SNR_sens_sqrt)/2;
        }
        
        //Se calcula el producto de los bias en log
        producto_bias_log=0;
        for (k=0;k<Nblancos;k++){
            if(indicador_blancos_pred[k+Nblancos*j]){
                producto_bias_log=producto_bias_log+log(Bias[j+Nparticulas*k]);
            }
        }
        
        pesos_pred[j]=Likelihood-producto_bias_log-log(pesos_existence_pred[j]); //Tambien hay que poner el logaritmo de pesos_existence_pred
        
    }
    
    
    //Hay que pasar a valor natural los pesos de las particulas
    peso_log_max=pesos_pred[0];
    for (j=0;j<Nparticulas;j++){
        if(pesos_pred[j]>peso_log_max){
            peso_log_max=pesos_pred[j];
        }
    }
    for (j=0;j<Nparticulas;j++){
        pesos_pred[j]=exp(-peso_log_max+pesos_pred[j]);
    }
    
    
    //Hay que normalizar los pesos de las particulas
    
    peso_total=0;
    for (j=0;j<Nparticulas;j++){
        peso_total=peso_total+pesos_pred[j];
    }
    for(j=0;j<Nparticulas;j++){
        pesos_pred[j]=pesos_pred[j]/peso_total;
    }
}
    
    
    
    
    
    
    
    
    
    
