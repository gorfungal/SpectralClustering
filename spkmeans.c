#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "spkmeans.h"

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

double** mat_mul(int n, double** m1, double** m2) {
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
                res[i][j] += m1[i][k] * m2[k][j];
            }
        }
    return res;
    }
}

void swap(double* xp, double* yp)
{
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
}
 
// A function to implement bubble sort
void bubbleSort(double arr[], int n)
{
    int i, j;
    for (i = 0; i < n - 1; i++){
        // Last i elements are already in place
        for (j = 0; j < n - i - 1; j++){
            if (arr[j] > arr[j + 1])
                swap(&arr[j], &arr[j + 1]);
        }
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


double** weighted_mat(double** allpoints, int n, int d){
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
				euclidean_norm+= pow(allpoints[i][k] - allpoints[j][k], (double)2);
			}
			euclidean_norm = math.sqrt(euclidean_norm);
			w_mat[i][j] = exp(-(euclidean_norm)/(double) 2));
		}
    }
    return w_mat;
}


double** diagonal_mat(int n, double** weighted){
    int i, j;
    double** w_mat = weighted;
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
    return diagonal_mat;
    /* need to free memory */
}


double** l_norm(int n, double** weighted){
    int i,j;
    double** I = (double**)malloc(n*sizeof(double*));
    assert(I);
    for (i=0; i<n; i++){
        I[i] = calloc(n * sizeof(double));
        assert(I[i]);
        I[i][i] = 1;
    }
    /* get diagonal ** -0.5 */
    double** diag = diagonal_mat(n, weighted);
    for (i=0; i<n; i++){
        diag[i][i] = 1.0 / pow(diag[i][i],0.5);
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
    return l_norm;
}

//                _   _               __       ___  
//               | | (_)             /_ |     |__ \ 
//  ___  ___  ___| |_ _  ___  _ __    | |        ) |
// / __|/ _ \/ __| __| |/ _ \| '_ \   | |       / / 
// \__ \  __/ (__| |_| | (_) | | | |  | |  _   / /_ 
// |___/\___|\___|\__|_|\___/|_| |_|  |_| (_) |____|
//                                                  
//                                                  

// where to calculate a_off_diag_sum?
void p_rot_mat(double** a_mat, double** v_mat, double* a_off_diag_sum, int* current_iter, int n){
    int i=0, j=0, i_max=0, j_max=0;
    double s=0, c=0, t=0, teta=0, sign_teta=0;
    double a_i_j=0, a_r_i=0, a_i_i=0, a_j_j=0, a_r_j=0;
    double a_i_j_max=0;
	double** p_mat=0;
    p_mat = (double **)malloc(n * sizeof(double *));
    if(!p_mat) printf("Error");
    for (i=0; i<n; i++){
        p_mat[i] = (double *)calloc(n, sizeof(double));	
        if(!p_mat[i]) printf("Error in p[%d]", i);
    }
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
	teta= (a_mat[j_max][j_max] - a_mat[i_max][i_max])/2*a_mat[i_max][j_max];
	if(teta >=0)
		sign_teta = 1;
	else
		sign_teta = -1;
	t = (sign_teta)/(fabs(teta) + sqrt(teta*teta+1));
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
        *a_off_diag_sum = new_a_off_diag_sum;
	}
    double** prev_a = a_mat;
	*a_mat = new_a_mat;
	for (i=0; i<n; i++){
        free(prev_a[i]);	
    }
	free(prev_a);
}

double** jacobi(double** lnorm, int n){
    int i, j;
    // init v_mat
    double** v_mat;
    v_mat = (double **)malloc(n * sizeof(double *));
    if(!v_mat) printf("Error");
    for (i=0; i<n; i++){
        v_mat[i] = (double *)calloc(n, sizeof(double));	
        if(!v_mat[i]) printf("Error in v[%d]", i);
    }
    double a_off = 0.0;
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            if(j!=i){
                a_off += pow(lnorm[i][j], 2.0);
            }
        }
    }

    int current_iter = 0;
    while(current_iter < MAX_ROT){
        p_rot_mat(lnorm, v_mat, &a_off, &current_iter, n);
        current_iter++;
    }
    // lnorm is changed in place
    return v_mat;
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

int Eigengap_heuristic(double** a_mat,double** v_mat, int n){
	double* eigenvalues = {0};
	double eigengap_max,eigengap_current = 0;
	for (i=0; i<n; i++){
		eigenvalues[i] = a_mat[i][i];
	}
	bubbleSort(eigenvalues);
	for (i=0; i < (int) ((n/2)+1); i++){ // 5/2 + 1 = 3 i will run on 1,2 as needed, 6/2 +1 = 4, i will run on 1,2,3 as needed. i=1,...,n/2
		double eigengap_current = abs(eigenvalues[i]-eigenvalues[i+1]);
		if(eigengap_current > eigengap_max){
			eigengap_max = eigengap_current;
		}
	}
	return (int) eigengap_max;
}
	
//  __  __       _       
// |  \/  |     (_)      
// | \  / | __ _ _ _ __  
// | |\/| |/ _` | | '_ \ 
// | |  | | (_| | | | | |
// |_|  |_|\__,_|_|_| |_|
//                      
double** readfile(char* filename, int* n, int* d){
    // read dimensions & number of points
    file = fopen(filename, "r");
    if (file == NULL){
        printf("An Error Has Occurred\n");
        return 1;
    }
    while (fscanf(file, "%lf", &num) != EOF){
        if ((c=getc(file)) == '\n'){
            *n++;
            *d++;
        }
        else if (c == ','){
            *d++;
        }
        else if (c == EOF){
            *n++;
            *d++;
            break;
        }
        else{
            *n++;
            *d++;
            break;
        }
    }
    *d /= *n;
    fseek(file, 0L, SEEK_SET);
    // init points array
    double** points = calloc(n, sizeof(double*));
    if(!points) printf("Problem with 'points' initialization");
    for(i=0; i<n; i++){
        points[i] = calloc(d, sizeof(double));
        if(!points[i]) printf("Problem with points%d", i);
    }
    // read points from file
    double num = 0.0; 
    int dim = 0; 
    int index = 0;
    while(fscanf(file, "%lf", &num) != EOF){
        if ((c=getc(file)) == '\n'){
            points[index][dim] = num;
            index++;
            dim=0;   
        }
        else if (c == ','){
            points[index][dim] = num;
            dim++;
        }
        else {
            points[index][dim] = num;
            break;
        }
    }
    fclose(file);
    return points;
}

void print_matrix(double** mat, int n, int k){
	for (i=0; i<n; i++){
        for (j=0; j<k; j++){
			if(j != k-1):
				printf("%0.4f ", mat[n][k]);
			else:
				printf("%0.4f\n", mat[n][k]);
		}
	}
}

double** spkmean(char* goal, char* filename, int k){
    int n=0; int d=0; int i, j;
    double** points = readfile(filename, &n, &d);
    double** w_mat;
    w_mat = weighted_mat(points, n, d);
    if(strcmp(goal, "wam") == 0){
        print_matrix(w_mat, n, n);
        return w_mat;
    }
    double** diag;
    diag = diagonal_mat(n, w_mat);
    if(strcmp(goal, "ddg") == 0){
        print_matrix(diag, n, n);
        return diag;
    }
    double** lnorm;
    lnorm = l_norm(n, w_mat);
    if(strcmp(goal, "lnorm") == 0){
        print_matrix(l_norm, n, n);
        return l_norm;
    }
    double** v_mat;
    // L_Norm -> A (Eigenvalues of L_Norm)
    v_mat = jacobi(lnorm, n);
    if(strcmp(goal, "jacobi") == 0){
        print_matrix(lnorm, n, n);
        print_matrix(v_mat, n, n);
        return lnorm;
    }
    // choosing k
    if(k == 0) k = Eigengap_heuristic(lnorm, v_mat, n);

    double** spk = (double**)malloc(k * sizeof(double*));
    if(!spk) printf("Error with spk");
    for(i=0; i<k; i++){
        spk[i] = v_mat[i];
    }
    // Init U and Transpose spk -> U
    double** U;
    U = (double **)malloc(n * sizeof(double *));
    if(!U) printf("Error with U");
    for (i=0; i<n; i++){
        U[i] = (double *)calloc(k, sizeof(double));	
        if(!U[i]) printf("Error in p[%d]", i);
    }
    for(i=0; i<k; i++){
        for(j=0; j<n; j++){
            U[j][i] = spk[i][j];
        }
    }
    // U -> T
    double row_sum;
    for(i=0; i<n; i++){
        row_sum = 0.0;
        for(j=0; j<k; j++){
            row_sum += pow(U[i][j], 2.0);
        }
        row_sum = pow(row_sum, 0.5);
        for(j=0; j<k; j++){
            U[i][j] = U[i][j] / row_sum;
        }
    }
    for(i=0; i<n; i++){
        free(w_mat[i]);
        free(diag[i]);
        free(lnorm[i]);
        free(v_mat[i]);
    }
    free(w_mat);
    free(diag);
    free(lnorm);
    free(v_mat);
    for(i=0; i<k; i++){
        free(spk[i]);
    }
    free(spk);
    return U;
}

int main(int argc, char *argv[]){
    int n = 0;
    int d = 0;
    int i = 0;
    char *goal = char[7];
    char *input_text = char[100];

    if(argc == 3){
        goal = argv[1];
        input_text = argv[2];
    }
    else printf("Invalid input!");
    if(strcmp(goal, "wam") == 0 || strcmp(goal, "ddg")==0 || strcmp(goal, "lnorm")==0 || strcmp(goal, "jacobi")==0){
        spkmean(goal, input_text);
    }

}