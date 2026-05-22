#include "mex.h"
#include "matrix.h"


#include <math.h>
//Subrutina para el problema de múltiple detección con etiquetas
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    //Definicion de variables de entrada
    double *particulas_act, *indicador_blancos_act,*Existence_matrix,*index_mult_existence;
    //Definición de variables intermedias
    int ndim_1, dims_1[2];
    int Nparticulas,Nblancos,Nexistence;
    int j,k;
    int index_existence,cardinality;
    double peso_part;
    
    //Definición de variables de salida
    double *prob_hypothesis,*prob_cardinality,*indice_particula_existence,*Npart_existence;
    
    //Inputs
    particulas_act=mxGetPr(prhs[0]);
    indicador_blancos_act=mxGetPr(prhs[1]);
    Existence_matrix=mxGetPr(prhs[2]);
    index_mult_existence=mxGetPr(prhs[3]);
    
    
    
    
    
    //Dimensiones vectores
    Nparticulas=mxGetN(prhs[0]);
    Nblancos=mxGetM(prhs[0])/4;
    Nexistence=mxGetN(prhs[2]);
    
    
    //Outputs
    ndim_1=2;
    dims_1[0] =1;
    dims_1[1] = Nexistence;
    plhs[0]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    prob_hypothesis = mxGetPr(plhs[0]);
    
    ndim_1=2;
    dims_1[0]=1;
    dims_1[1]=Nblancos+1;
    plhs[1]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    prob_cardinality=mxGetPr(plhs[1]);
    
     ndim_1=2;
    dims_1[0]=1;
    dims_1[1]=Nparticulas;
    plhs[2]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    indice_particula_existence=mxGetPr(plhs[2]);
    
    ndim_1=2;
    dims_1[0]=1;
    dims_1[1]=Nexistence;
    plhs[3]=mxCreateNumericArray(ndim_1, dims_1, mxDOUBLE_CLASS, mxREAL);
    Npart_existence=mxGetPr(plhs[3]);
    
    
    
    //Reserva de memoria
    
    
    
    //Operaciones
    peso_part=(double)1/Nparticulas;
    
    
    for (j=0;j<Nexistence;j++){
        prob_hypothesis[j]=0;
        Npart_existence[j]=0;
    }
    
    for (j=0;j<Nblancos+1;j++){
        prob_cardinality[j]=0;
    }
    
    
    
    for (j=0;j<Nparticulas;j++){
        
        index_existence=0;
        cardinality=0;
        for (k=0;k<Nblancos;k++){
            
            index_existence=index_existence+indicador_blancos_act[k+Nblancos*j]*index_mult_existence[k];
            
            if(indicador_blancos_act[k+Nblancos*j]==1){
                cardinality=cardinality+1;
            }
        }
        
        indice_particula_existence[j]=index_existence;
        prob_hypothesis[index_existence]=prob_hypothesis[index_existence]+peso_part;
        Npart_existence[index_existence]=Npart_existence[index_existence]+1;
        prob_cardinality[cardinality]=prob_cardinality[cardinality]+peso_part;
        
    }
    
    
    
}
