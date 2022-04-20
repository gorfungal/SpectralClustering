#include <math.h>
#include <assert.h>
#include <stdlib.h>

EPS = 10**-15;
MAX_ROT = 100;




//   _____                           _ 
//  / ____|                         | |
// | |  __  ___ _ __   ___ _ __ __ _| |
// | | |_ |/ _ \ '_ \ / _ \ '__/ _` | |
// | |__| |  __/ | | |  __/ | | (_| | |
//  \_____|\___|_| |_|\___|_|  \__,_|_|
//                                     
//                                     

double*** mat_mul(int n, double*** m1, double*** m2) {
    /* don't think memory freeing is needed */
    int i,j;
    double** res = (double**)malloc(n*sizeof(double*));
    assert(res);
    for (i=0; i<n; i++){
        res[i] = (double*)malloc(n*sizeof(double));
        assert(res[i]);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            res[i][j] = 0.0;
            for (int k = 0; k < n; k++) {
                res[i][j] += (*m1)[i][k] * (*m2)[k][j];
            }
        }
    return &res;
    }
}



//                _   _               __       __ 
//               | | (_)             /_ |     /_ |
//  ___  ___  ___| |_ _  ___  _ __    | |      | |
// / __|/ _ \/ __| __| |/ _ \| '_ \   | |      | |
// \__ \  __/ (__| |_| | (_) | | | |  | |  _   | |
// |___/\___|\___|\__|_|\___/|_| |_|  |_| (_)  |_|
//                                                
//                                                


double*** weighted_mat(double** allpoints, int n, int d){
    int i, j, k;
	double euclidean_norm;
    double** w_mat = (double **)malloc(n * sizeof(double *));
    assert(w_mat);
    for (i = 0; i<n; i++){
        w_mat[i] = (double *)calloc(n, sizeof(double));
        assert(w_mat[i]);
    }
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            euclidean_norm = 0;
			for (k=0; k<d; k++){
				euclidean_norm+= pow((double) 2,allpoints[i][k] - allpoints[j][k]);
			}
			euclidean_norm = math.sqrt(euclidean_norm)
			w_mat[i][j] = exp(-(euclidean_norm)/(double) 2));
		}
    }
    return &w_mat;
}


double*** diagonal_mat(int n, double*** weighted){
    int i, j;
    double** w_mat = * weighted;
    double** diag = (double **)malloc(n * sizeof(double *));
    assert(diag);
    for (i = 0; i<n; i++){
        diag[i] = (double *)calloc(n, sizeof(double));
        assert(diag[i]);
    }
    for (i=0; i<n; i++){
        double sum = 0.0;
        for (j=0; j<n; j++){
            sum += w_mat[j][i];
        }
        diag[i][i] = sum;
    }
    return &diagonal_mat;
}

double quake_exp(double num){
    /* does double work the same as float? */
    long i;
    double x2, y;
    const double t_h = 1.5;
    
    x2 = num * 0.5;
    y = num;
    i = * (long*)&y;
    i = 0x5f3759df - (i>>1);
    y = y*(t_h - (x2*y*y));
    y = y*(t_h - (x2*y*y)); /* necessary? */
    return y;
}

double*** l_norm(int n, double *** weighted){
    int i,j;
    double** I = (double**)malloc(n*sizeof(double*));
    assert(I);
    for (i=0; i<n; i++){
        I[i] = calloc(n * sizeof(double));
        assert(I[i]);
        I[i][i] = 1;
    }
    /* get diagonal ** -0.5 */
    double** diag = *diagonal_mat(n, weighted);
    for (i=0; i<n; i++){
        /* check precision */
        diag[i][i] = quake_exp(diag[i][i]);
    }

    /* D^-0.5 * W * D^-0.5 */
    double** dwd = *mat_mul(n, &diag, weighted);
    dwd = *mat_mul(n, &dwd, &diag);

    double** l_norm = (double **)malloc(n*sizeof(double*));
    assert(l_norm);
    for (i=0; i<n; i++){
        l_norm[i] = calloc(n * sizeof(double));
        assert(l_norm[i]);
        for (j=0; j<n; j++){
            l_norm[i][j] = I[i][j] - dwd[i][j];
        }
    }
    return &l_norm;
}

//                _   _               __       ___  
//               | | (_)             /_ |     |__ \ 
//  ___  ___  ___| |_ _  ___  _ __    | |        ) |
// / __|/ _ \/ __| __| |/ _ \| '_ \   | |       / / 
// \__ \  __/ (__| |_| | (_) | | | |  | |  _   / /_ 
// |___/\___|\___|\__|_|\___/|_| |_|  |_| (_) |____|
//                                                  
//                                                  

void p_rot_mat(double** a_mat, double** v_mat, double* a_off_diag_sum, int* current_iter, int n){
    int i=0, j=0, i_max=0, j_max=0;
    double s=0, c=0, t=0, teta=0, sign_teta=0;
    double a_i_j=0, a_r_i=0, a_i_i=0, a_j_j=0, a_r_j=0;
    double a_i_j_max=0;
	double** p_mat=0;
	//find max
	for (i=0; i<n; i++){
        for (j=i+1; j<n; j++){
			if(fabs(a_mat[i][j]) > fabs(a_i_j_max)){
				a_i_j_max =a_mat[i][j]; 
				i_max = i; 
				j_max = j;
			}
		}
    }
	//calculate variables for p mat
	teta = (a_mat[j_max][j_max] - a_mat[i_max][i_max])/2*a_mat[i_max][j_max];
    // why is it 1 for both cases?
	if(teta >=0)
		sign_teta = 1;
	else
		sign_teta = 1;
	t = (sign_teta)/fabs(teta) + sqrt(teta*teta+1);
	c = 1/(t*t + 1);
	s= t*c;
	//p mat calculation
	for (i=0; i<n; i++){
        for (j=0; j<n; j++){
			if(i==j)
				p_mat[i][j] = 1; 
			if(i == i_max && j==j_max)
				p_mat[i][j] = s; 
			if(i == j_max && j==j_max)
				p_mat = c; 
			if(i == i_max && j==i_max)
				p_mat = c;
			if(j == i_max && i==j_max)
				p_mat = -s;
		}
	}
	//v mat calculation
	v_mat = mat_mul(n, v_mat, p);
	//a_mat calculation
    double** new_a_mat = (double **)malloc(n * sizeof(double *));
	double new_a_off_diag_sum=0;
	for (i=0; i<n; i++){
        new_a_mat[i] = (double *)calloc(n, sizeof(double));	
        double temp_a_value=0;
		for (j=i+1; j<n; j++){ // check if correct
			if(i == i_max && j == j_max){
				temp_a_value = 0;
				new_a_mat[i][j] = temp_a_value;
				new_a_mat[i][j] = temp_a_value;
			}
			else if(i == i_max && j==i_max){
				temp_a_value = c*c*a_mat[i][i] + s*s*a[j_max][j_max] -2*s*c*a_mat[i][j_max];
				new_a_mat[i][j] = temp_a_value; 
				new_a_mat[j][i] = temp_a_value;
			}
			else if(i == j_max && j==j_max){
				temp_a_value = c*c*a_mat[j][j] + s*s*a[i_max][i_max] +2*s*c*a_mat[i_max][j];
				new_a_mat[i][j] = temp_a_value; 
				new_a_mat[j][i] = temp_a_value;
			}
			else if(j == i_max){
				temp_a_value = c*a_mat[j][i_max] - s*a_mat[j][j_max];
				new_a_mat[i][j] = temp_a_value; 
				new_a_mat[j][i] = temp_a_value;
			}
			else if(j == i_max){
				temp_a_value = c*a_mat[j][j_max] + s*a_mat[j][i_max];
				new_a_mat[i][j] = temp_a_value; 
				new_a_mat[j][i] = temp_a_value;
            }
			else {
				temp_a_value = a_mat[i][j];
				new_a_mat[i][j] = temp_a_value; 
				new_a_mat[j][i] = temp_a_value; 
			}
			new_a_off_diag_sum += new_a_mat[j][i]*new_a_mat[j][i] + new_a_mat[i][j]*new_a_mat[i][j]; 	
		}
	}
	if(*a_off_diag_sum-new_a_off_diag_sum <= EPS){
		*current_iter = MAX_ROT;
	}
	*a_mat = new_a_mat;
	for (i=0; i<n; i++){
        free(new_a_mat[i]);	
    }
	free(new_a_mat);
}

//
//                _   _               __       ____  
//               | | (_)             /_ |     |___ \ 
//  ___  ___  ___| |_ _  ___  _ __    | |       __) |
// / __|/ _ \/ __| __| |/ _ \| '_ \   | |      |__ < 
// \__ \  __/ (__| |_| | (_) | | | |  | |  _   ___) |
// |___/\___|\___|\__|_|\___/|_| |_|  |_| (_) |____/ 
//                                                   
//
int compare(float num1, float num2){
    if(num1<num2) return -1;
    else if(num1==num2) return 0;
    else return 1;
}

int num_of_clusters(float **eigenvalues){
    int n = sizeof(*eigenvalues)/sizeof(float);
    qsort(*eigenvalues, n, sizeof(float), compare);
    int cnt = 0;
    int argmax = 0;
    float delta;
    float delta_max = 0;
    for(int i=1; i < (n/2); i++){
        delta = *eigenvalues[i] - *eigenvalues[i-1];
        if(delta > delta_max){
            delta_max = delta;
            argmax = i;
        }
    }
    return argmax;
}


//  __  __       _       
// |  \/  |     (_)      
// | \  / | __ _ _ _ __  
// | |\/| |/ _` | | '_ \ 
// | |  | | (_| | | | | |
// |_|  |_|\__,_|_|_| |_|
//                      
