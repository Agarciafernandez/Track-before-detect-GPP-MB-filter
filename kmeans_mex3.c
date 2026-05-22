#include "mex.h"
#include "matrix.h"
#include "matrixutils.h"
#include "statutils.h"
#include <math.h>

//Funcion obtenida de la carpeta Tracking Nacimiento muerte. Hago alguna ligera modificacion
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //Definicion de variables de entrada
    double *particulas_pred, *pesos_pred,*Perm_matrix,*Nmax_iteraciones;
    //Definición de variables intermedias
    int flag, i,ndim_1,dims_1[2],k,q;
    int Nparticulas,Nblancos,Nperm;
    double resta_pos_x,resta_vel_x,resta_pos_y,resta_vel_y;
    int q_minimo;
    double coste_minimo,coste_q;
    double *particula_aux;
    int contador,index_perm;
    
    //Definición de variables de salida
    double *particulas_sal, *Estado_med_blancos, *indicador_cambio;
    //Inputs
    particulas_pred=mxGetPr(prhs[0]);
    pesos_pred=mxGetPr(prhs[1]);
    Perm_matrix=mxGetPr(prhs[2]);
    Nmax_iteraciones=mxGetPr(prhs[3]);
    
    //Dimensiones vectores
    Nparticulas=mxGetN(prhs[0]);
    Nblancos=mxGetM(prhs[0])/4;
    Nperm=mxGetM(prhs[2]);
    
    
    //Outputs
    ndim_1=2;
    dims_1[0] = 4*Nblancos;
    dims_1[1] =  Nparticulas;
    plhs[0]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    particulas_sal = mxGetPr(plhs[0]);
    ndim_1=2;
    dims_1[0] = 4*Nblancos;
    dims_1[1] = 1;
    plhs[1]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    Estado_med_blancos = mxGetPr(plhs[1]);
    
    ndim_1=1;
    dims_1[0] = 1;
    
    plhs[2]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    indicador_cambio = mxGetPr(plhs[2]);
    
    //Reserva memoria
    particula_aux=(double *) mxCalloc(4*Nblancos, sizeof(double));
    
    //Operaciones
    flag=0;//Cuando sea 1 el algoritmo ha terminado
    for(i=0;i<4*Nblancos*Nparticulas;i++){
        particulas_sal[i]=particulas_pred[i];
    }
    
    contador=0;
    indicador_cambio[0]=0;
    while (flag==0 && contador<Nmax_iteraciones[0]){
        flag=1;
        
        for (i=0;i<4*Nblancos;i++){
            Estado_med_blancos[i]=0;
        }
        
        for (i=0;i<Nblancos;i++){
            
            for (k=0;k<Nparticulas;k++){
                Estado_med_blancos[4*i]=Estado_med_blancos[4*i]+pesos_pred[k]*particulas_sal[4*i+4*Nblancos*k];
                Estado_med_blancos[4*i+1]=Estado_med_blancos[4*i+1]+pesos_pred[k]*particulas_sal[4*i+4*Nblancos*k+1];
                Estado_med_blancos[4*i+2]=Estado_med_blancos[4*i+2]+pesos_pred[k]*particulas_sal[4*i+4*Nblancos*k+2];
                Estado_med_blancos[4*i+3]=Estado_med_blancos[4*i+3]+pesos_pred[k]*particulas_sal[4*i+4*Nblancos*k+3];
                
            }
        }
        for (k=0;k<Nparticulas;k++){
            coste_q=0;
            for (i=0;i<Nblancos;i++){
                resta_pos_x=Estado_med_blancos[4*i]-particulas_sal[4*i+4*Nblancos*k];
                resta_vel_x=Estado_med_blancos[4*i+1]-particulas_sal[4*i+4*Nblancos*k+1];
                resta_pos_y=Estado_med_blancos[4*i+2]-particulas_sal[4*i+4*Nblancos*k+2];
                resta_vel_y=Estado_med_blancos[4*i+3]-particulas_sal[4*i+4*Nblancos*k+3];
                coste_q=coste_q+resta_pos_x*resta_pos_x+resta_vel_x*resta_vel_x+resta_pos_y*resta_pos_y+resta_vel_y*resta_vel_y;               
            }
            
            coste_minimo=coste_q;
            q_minimo=0;
            for (q=1;q<Nperm;q++){
                coste_q=0;
                for (i=0;i<Nblancos;i++){
                    index_perm=(int)Perm_matrix[q+Nperm*i];
                    resta_pos_x=Estado_med_blancos[4*i]-particulas_sal[4*index_perm+4*Nblancos*k];
                    resta_vel_x=Estado_med_blancos[4*i+1]-particulas_sal[4*index_perm+4*Nblancos*k+1];
                    resta_pos_y=Estado_med_blancos[4*i+2]-particulas_sal[4*index_perm+4*Nblancos*k+2];
                    resta_vel_y=Estado_med_blancos[4*i+3]-particulas_sal[4*index_perm+4*Nblancos*k+3];
                    coste_q=coste_q+resta_pos_x*resta_pos_x+resta_vel_x*resta_vel_x+resta_pos_y*resta_pos_y+resta_vel_y*resta_vel_y;
                    
                }
                if(coste_q<coste_minimo){
                    coste_minimo=coste_q;
                    q_minimo=q;
                }
                
            }
            if(q_minimo>0){
                //Ha habido cambio->Se pone el flag a cero y hay que hacer otra permutacion
                flag=0;
                indicador_cambio[0]=1;

                
                
                for (i=0;i<Nblancos;i++){
                    particula_aux[4*i]=particulas_sal[4*i+4*Nblancos*k];
                    particula_aux[4*i+1]=particulas_sal[4*i+4*Nblancos*k+1];
                    particula_aux[4*i+2]=particulas_sal[4*i+4*Nblancos*k+2];
                    particula_aux[4*i+3]=particulas_sal[4*i+4*Nblancos*k+3];
                }
                
                
                for (i=0;i<Nblancos;i++){
                    index_perm=(int)Perm_matrix[q_minimo+Nperm*i];
                    //Se permuta el indicador_blancos_sal y la particulas
                    particulas_sal[4*index_perm+4*Nblancos*k]=particula_aux[4*i];
                    particulas_sal[4*index_perm+4*Nblancos*k+1]=particula_aux[4*i+1];
                    particulas_sal[4*index_perm+4*Nblancos*k+2]=particula_aux[4*i+2];
                    particulas_sal[4*index_perm+4*Nblancos*k+3]=particula_aux[4*i+3];
                }
            }
        }
        contador=contador+1;
    }
    //printf("%d\n",contador);
}



