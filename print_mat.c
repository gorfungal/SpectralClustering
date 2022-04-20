void print_mat(double** mat,int n,int k){
	for (i=0; i<n; i++){
        for (j=0; j<k; j++){
			if(j != k-1):
				printf("%0.4f ", mat[n][k]);
			else:
				printf("%0.4f\n", mat[n][k]);
		}
	}
}