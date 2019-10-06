
void classifier_MMM2(int (*p1), int (*K1), int (*n1), double *x1, double *alpha, double *Pi1, double *gamma1, int *id, int (*missingStates), int (*lenMissing1)){

	int i, j, r, p, K, n, lenMissing; 
	double ***x; //, **nj;
	
	double **gamma, ***MM, ***MM_rowOrder; 
	
	p = (*p1);
	K = (*K1);
	n = (*n1);
	lenMissing = (*lenMissing1);

	MAKE_3ARRAY(x, n, p, p);	
	MAKE_3ARRAY(MM, p, p, K);
	MAKE_MATRIX(gamma, n, K);
	MAKE_3ARRAY(MM_rowOrder, K, p, p);
	
	//Rprintf("\nTotal missing states : %d", lenMissing);
	//Rprintf("\nVal p : %d", p);

	if(lenMissing > 0)
		array1to3_missState(n, p, p, x1, x, missingStates, lenMissing);	 // consider missing state indexes while filling the x array from vector
	else 
		array1to3_(n, p, p, x1, x);
	
	array1to3(p, p, K, Pi1, MM);
	array1to3_(K, p, p, Pi1, MM_rowOrder);
	array1to2(n, K, gamma1, gamma);	
	
	

/*	Rprintf("\n\n values x ... classifier_MMM2 (libClickClust.c) \n");
	
	for (i=0; i<n; i++){
	  for (j=0; j<p; j++){
	    for (r=0; r<p; r++){
	      Rprintf("%f ", x[i][j][r]);
	    }
	    Rprintf("\n");
	  }
	  Rprintf("\n");
	}
*/

	classify2(p, n, x, alpha, MM, MM_rowOrder, gamma, id, K);
	

	array3to1(p, p, K, Pi1, MM);
	array3to1_(K, p, p, Pi1, MM_rowOrder);
	array2to1(n, K, gamma1, gamma);
	
	
	FREE_3ARRAY(MM);
	FREE_3ARRAY(MM_rowOrder);
	FREE_MATRIX(gamma);	

}
