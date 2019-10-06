
#include <R.h>
#include<math.h>
#include "array.h"
#define inf 1e+40;
#include "ClickClust.h"

void classify2(int p, int n, double ***x, double *alpha, double ***Pi, double ***Pi_rowOrder, double **gamma, int *id, int K){
	
	int k, i, j, r;
	double gmax;
	
	// assign equal weight to each mixture
	for (k=0; k<K; k++){
			alpha[k] = 1.0 / K;
	}
		
	//Rprintf("\n\n values Pi_rowOrder ... classify (libEM.c) \n");
	
/*	for (k=0; k<K; k++){
	  for (j=0; j<p; j++){
	    for (r=0; r<p; r++){
	      Rprintf("%f ", Pi_rowOrder[k][j][r]);
	    }
	    Rprintf("\n");
	  }
	  Rprintf("\n");
	}

	Rprintf("\n\n values x ... classify (libEM.c) \n");
	
	for (i=0; i<n; i++){
	  for (j=0; j<p; j++){
	    for (r=0; r<p; r++){
	  	  //x[i][j][r] = (long double)(x[i][j][r]);
	      Rprintf("%f ", x[i][j][r]);
	    }
	    Rprintf("\n");
	  }
	  Rprintf("\n");
	}
*/

	Estep2(p, n, x, alpha, Pi_rowOrder, gamma, K);

	for (i=0; i<n; i++){
		gmax = gamma[i][0];
		id[i] = 0;
		for (k=1; k<K; k++){
			if (gamma[i][k] > gmax){
				gmax = gamma[i][k];
				id[i] = k;
			}
		}
		//Rprintf("\n\nAssigned cluster id by posterior_ik to seq. %d = %d", i, id[i]);
	}
}


void Estep2(int p, int n, double ***x, double *alpha, double ***Pi_rowOrder, double **gamma, int K){
	
	// computing posterior of sequene i with each cluster separately and divide each value by denominator, totalDenom[i] 

	long double **posterior_ik; 			
	long double *totalDenom, prod1;
	double prodTemp, *totalDenomTemp; 		// store results in gamma ... in double

	MAKE_MATRIX(posterior_ik, n, K);		// make a new matrix exactly as gamma
	MAKE_1ARRAY(totalDenom,n); 				// to store denomenator for n sequences in long double
	MAKE_1ARRAY(totalDenomTemp,n); 			// to store denomenator for n sequences in double
	
	int i, k, j, r;  

	// computing denom... likelihood for all clusters
	for(i=0; i<n; i++){
		totalDenom[i] = 0.0;
		totalDenomTemp[i] = 0.0;
		
		for (k=0; k<K; k++){
			posterior_ik[i][k] = 0.0;
			gamma[i][k] = 0.0; 
			prod1 = 1.0;
			prodTemp = 0.0; 
			for (j=0; j<p; j++){
			    for (r=0; r<p; r++){

					prod1 = prod1 * powl((long double)(Pi_rowOrder[k][j][r]),(long double)(x[i][j][r])); 
	    			prodTemp = prodTemp + x[i][j][r] * log(Pi_rowOrder[k][j][r]);
			    }
		  	}
		  	prod1 = (long double) (alpha[k]) * prod1; 	
			prodTemp = exp(prodTemp);
			prodTemp = alpha[k] * prodTemp; 	
			posterior_ik[i][k] = prod1; 
			totalDenom[i] = totalDenom[i] + posterior_ik[i][k];

			gamma[i][k] = prodTemp; 
			totalDenomTemp[i] = totalDenomTemp[i] + gamma[i][k];
		}
		//Rprintf("\n=========================");	
		//Rprintf("\n\n totalDenom for seq.i (%d) = %Lf ", i,totalDenom[i]);
	}


	for(i=0; i<n; i++){
		for (k=0; k<K; k++){
			posterior_ik[i][k] = posterior_ik[i][k] / totalDenom[i];

			gamma[i][k] = gamma[i][k] / totalDenomTemp[i];
		}
	}

	FREE_MATRIX(posterior_ik);	
	FREE_1ARRAY(totalDenom); 				
	FREE_1ARRAY(totalDenomTemp); 				
}

