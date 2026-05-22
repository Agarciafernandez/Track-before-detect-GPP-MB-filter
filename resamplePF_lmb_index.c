#include "mex.h"
#include "matrix.h"
#include "matrixutils.h"
#include "statutils.h"
#include <math.h>

// Calcula el resample de las partículas
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    
    //Definicion variables de entrada
    double *particulas_pred, *pesos_pred, *index_parent_pred,*indicador_blancos_pred;
    //Definicion variables intermedias
    int dims_1[2];
    double *cdf;
    double u1, uj;
    int k, ir, jr, Nparticulas, Ndim,Nblancos, i;
    //Definicion variables de salida
    double *particulas_act, *pesos_act, *index_parent_act,*indicador_blancos_act;
    
    //Inputs
    particulas_pred=mxGetPr(prhs[0]);
    pesos_pred=mxGetPr(prhs[1]);
    index_parent_pred=mxGetPr(prhs[2]);
    indicador_blancos_pred=mxGetPr(prhs[3]);
    
    
    //Dimensiones vectores
    
    Nparticulas=mxGetN(prhs[0]);
    Ndim=mxGetM(prhs[0]);
    Nblancos=Ndim/4;
    
    
    //Outputs
    dims_1[0]=Ndim;
    dims_1[1]=Nparticulas;
    plhs[0]=mxCreateNumericArray(2, dims_1, mxDOUBLE_CLASS, mxREAL);
    particulas_act = mxGetPr(plhs[0]);
    
    dims_1[0]=Nparticulas;
    dims_1[1]=Nblancos;
    plhs[1]=mxCreateNumericArray(2, dims_1, mxDOUBLE_CLASS, mxREAL);
    index_parent_act = mxGetPr(plhs[1]);
    
    dims_1[0]=Nblancos;
    dims_1[1]=Nparticulas;
    plhs[2]=mxCreateNumericArray(2, dims_1, mxDOUBLE_CLASS, mxREAL);
    indicador_blancos_act = mxGetPr(plhs[2]);
    
    
    
    //Reserva de memoria
    cdf=(double *)mxCalloc(Nparticulas, sizeof(double));
    
    
    //Se calcula la cdf
    
    cdf[0]=pesos_pred[0];
    
    //Este bucle empieza en 1
    for (k=1;k<Nparticulas;k++){
        cdf[k]=cdf[k-1]+pesos_pred[k];
        
    }
    u1=unirnd(0, 1)/Nparticulas;
    ir=0;
    for (jr=0;jr<Nparticulas;jr++){
        uj=u1+((double)jr)/Nparticulas;
        
        while(uj>cdf[ir]){
            ir=ir+1;
        }
        for(i=0;i<Ndim;i++){
            particulas_act[i+Ndim*jr]=particulas_pred[i+Ndim*ir];
        }
        for (i=0;i<Nblancos;i++){
            index_parent_act[jr+Nparticulas*i]=index_parent_pred[ir+Nparticulas*i];           
            indicador_blancos_act[i+Nblancos*jr]=indicador_blancos_pred[i+Nblancos*ir];
        }
        
     
    }
}