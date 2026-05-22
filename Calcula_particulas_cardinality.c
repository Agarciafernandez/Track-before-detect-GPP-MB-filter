#include "mex.h"
#include "matrix.h"


#include <math.h>
//Subrutina para el problema de múltiple detección con etiquetas
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //Definicion de variables de entrada
    double *particulas_act, *indicador_blancos_act,*cardinality_est,*indices_card,*Nparticulas_card;
    //Definición de variables intermedias
    int ndim_1, dims_1[2];
    int Nparticulas,Nblancos,Nblancos_sal;
    int j,k;
    int index_salida,index_blanco_salida;
    
    //Definición de variables de salida
    double *particulas_act_card;
    
    //Inputs
    particulas_act=mxGetPr(prhs[0]);
    indicador_blancos_act=mxGetPr(prhs[1]);
    cardinality_est=mxGetPr(prhs[2]);
    indices_card=mxGetPr(prhs[3]);
    Nparticulas_card=mxGetPr(prhs[4]);
    
    
    
    //Dimensiones vectores
    Nparticulas=mxGetN(prhs[0]);
    Nblancos=mxGetM(prhs[0])/4;
    
    Nblancos_sal=(int)cardinality_est[0];
    
    //Outputs
    ndim_1=2;
    dims_1[0] = 4*Nblancos_sal;
    dims_1[1] = (int) Nparticulas_card[0];
    plhs[0]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    particulas_act_card = mxGetPr(plhs[0]);
    

    
    //Reserva de memoria
    

    //Operaciones
    index_salida=0;
    
    for (j=0;j<Nparticulas;j++){
        if(indices_card[j]==1){
            index_blanco_salida=0;
            
            for (k=0;k<Nblancos;k++){
                if(indicador_blancos_act[k+Nblancos*j]==1){
                    particulas_act_card[4*index_blanco_salida+4*Nblancos_sal*index_salida]=particulas_act[4*k+4*Nblancos*j];
                    particulas_act_card[4*index_blanco_salida+1+4*Nblancos_sal*index_salida]=particulas_act[4*k+1+4*Nblancos*j];
                    particulas_act_card[4*index_blanco_salida+2+4*Nblancos_sal*index_salida]=particulas_act[4*k+2+4*Nblancos*j];
                    particulas_act_card[4*index_blanco_salida+3+4*Nblancos_sal*index_salida]=particulas_act[4*k+3+4*Nblancos*j];
                    
                    index_blanco_salida=index_blanco_salida+1;                                 
                }
            }           
            index_salida=index_salida+1;            
        }      
    }       
}
