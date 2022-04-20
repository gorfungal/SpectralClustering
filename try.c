#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define MAX_ROT 100

//   _____                           _ 
//  / ____|                         | |
// | |  __  ___ _ __   ___ _ __ __ _| |
// | | |_ |/ _ \ '_ \ / _ \ '__/ _` | |
// | |__| |  __/ | | |  __/ | | (_| | |
//  \_____|\___|_| |_|\___|_|  \__,_|_|
//                                     
//                                     

double** matrix_allocation(n, k)
{
    double** mat;
    mat = malloc(n*sizeof(double*));
    if (!mat) printf("problem with memory allocation mat");
    for (int point_i = 0; point_i < n; point_i++) {
        mat[point_i] = calloc(k, sizeof(double));
    }
    return mat;
}

void free_matrix_allocation(double** mat, int n) {
    for (int point_i = 0; point_i < n; point_i++) {
        if(mat[point_i] != NULL){
            free(mat[point_i]);
        }
    }
    if (mat != NULL) {
        free(mat);
    }
}


void print_matrix(double** mat, int n, int k) {
    int i = 0, j = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < k; j++) {
            if (j != k - 1) {
                printf("%0.4f ", mat[i][j]);
            }
            else {
                printf("%0.4f\n", mat[i][j]);
            }
        }
    }
}

double** mat_mul(int n, double** m1, double** m2) {
    /* don't think memory freeing is needed */
    int i, j;
    double** res = matrix_allocation(n, n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            res[i][j] = 0.0;
            for (int k = 0; k < n; k++) {
                res[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return res;
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
    for (i = 0; i < n - 1; i++) {
        // Last i elements are already in place
        for (j = 0; j < n - i - 1; j++) {
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


double** weighted_mat(double** allpoints, int n, int d) {
    int i, j, k;
    double euclidean_norm=0,w_i_j=0;
    double** w_mat = matrix_allocation(n, n);
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            if (i == j) {
                w_mat[i][j] = 0;
                continue;
            }
            euclidean_norm = 0;
            for (k = 0; k < d; k++) {
                euclidean_norm += pow(allpoints[i][k] - allpoints[j][k], (double) 2);
            }
            euclidean_norm = sqrt(euclidean_norm);
            w_i_j = exp(-(euclidean_norm) / 2.0);
            if (w_i_j > 0) {
                w_mat[i][j] = w_i_j;
                w_mat[j][i] = w_i_j;
            }
            else{
                w_mat[i][j] = 0;
                w_mat[j][i] = 0;
            }
        }
    }
    return w_mat;
}


double** diagonal_mat(int n, double** weighted) {
    int i, j;
    double** w_mat = weighted;
    double** diag = matrix_allocation(n, n);
    for (i = 0; i < n; i++) {
        double sum = 0.0;
        for (j = 0; j < n; j++) {
            sum += w_mat[i][j];
        }
        diag[i][i] = sum;
        diag[i][i] = 1.0 / pow(diag[i][i], 0.5);
    }
    return diag;
}

double** l_norm(int n, double** weighted, double** diag) {
    int i, j;
    double** I = matrix_allocation(n, n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                I[i][i] = 1.0;
            }
            else {
                I[i][j] = 0.0;
            }
        }
    }
    /* I - D^-0.5 * W * D^-0.5 */
    double** l_norm = mat_mul(n, diag, weighted);
    l_norm = mat_mul(n, l_norm, diag);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            l_norm[i][j] = I[i][j] - l_norm[i][j];
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

void p_rot_mat(double** a_mat, double** v_mat, double* a_off_diag_sum, int* current_iter, int n) {
    int i = 0, j = 0, i_max = 0, j_max = 0;
    double s = 0, c = 0, t = 0, teta = 0, sign_teta = 0;
    double a_i_j = 0, a_r_i = 0, a_i_i = 0, a_j_j = 0, a_r_j = 0;
    double a_i_j_max = 0.0;
    double** p_mat = matrix_allocation(n, n);
    //unit matrix
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                p_mat[i][i] = 1.0;
            }
            else {
                p_mat[i][j] = 0.0;
            }
        }
    }
    //find max
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            if (fabs(a_mat[i][j]) > fabs(a_i_j_max)) {
                a_i_j_max = a_mat[i][j];
                i_max = i;
                j_max = j;
            }
        }
    }
    //calculate variables for p mat
    teta = (a_mat[j_max][j_max] - a_mat[i_max][i_max]) / 2 * a_mat[i_max][j_max];
    if (teta >= 0) {
        sign_teta = 1;
    }
    else {
        sign_teta = -1;
    }
    t = (sign_teta) / (fabs(teta) + sqrt((teta * teta) + 1));
    c = 1 / sqrt(((t * t) + 1));
    s = t * c;
    //p mat calculation
    printf("a_mat\n");
    print_matrix(a_mat, n, n);
    p_mat[i_max][j_max] = s;
    p_mat[i_max][i_max] = c;
    p_mat[j_max][j_max] = c;
    p_mat[j_max][i_max] = -s;
    printf("p_mat\n");
    print_matrix(p_mat, n, n);
    //v mat calculation
    double** new_v_mat = mat_mul(n, v_mat, p_mat);
    //a_mat calculation
    double** new_a_mat = matrix_allocation(n, n);
    double new_a_off_diag_sum = 0;
    for (i = 0; i < n; i++) {
        double temp_a_value = 0;
        for (j = i; j < n; j++) { // check if correct
            if (i == i_max && j == j_max) {
                temp_a_value = 0;
                new_a_mat[i][j] = temp_a_value;
                new_a_mat[j][i] = temp_a_value;
            }
            else if (i == i_max && j == i_max) {
                temp_a_value = c * c * a_mat[i_max][i_max] + s * s * a_mat[j_max][j_max] - 2 * s * c * a_mat[i_max][j_max];
                new_a_mat[i][j] = temp_a_value;
            }
            else if (i == j_max && j == j_max) {
                temp_a_value = c * c * a_mat[j_max][j_max] + s * s * a_mat[i_max][i_max] + 2 * s * c * a_mat[i_max][j];
                new_a_mat[i][j] = temp_a_value;
            }
            else if (j == i_max) {
                temp_a_value = c * a_mat[i][i_max] - s * a_mat[i][j_max];
                new_a_mat[i][j] = temp_a_value;
                new_a_mat[j][i] = temp_a_value;
            }
            else if (j == j_max) {
                temp_a_value = c * a_mat[i][j_max] + s * a_mat[i][i_max];
                new_a_mat[i][j] = temp_a_value;
                new_a_mat[j][i] = temp_a_value;
            }
            else {
                temp_a_value = a_mat[i][j];
                new_a_mat[i][j] = temp_a_value;
                new_a_mat[j][i] = temp_a_value;
            }
            new_a_off_diag_sum += 2 * pow(new_a_mat[j][i], 2.0);
        }
    }
    printf("new_a_mat\n");
    print_matrix(new_a_mat, n, n);
    if (*a_off_diag_sum - new_a_off_diag_sum <= pow(10.0, -15.0)) {
        *current_iter = MAX_ROT;
        
    }
    a_off_diag_sum = &new_a_off_diag_sum;
    //copy to a mat
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) { // check if correct
            a_mat[i][j] = new_a_mat[i][j];
        }
    }
    //copy to a mat
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) { // check if correct
            v_mat[i][j] = new_v_mat[i][j];
        }
    }
    free_matrix_allocation(new_v_mat, n);
    free_matrix_allocation(new_a_mat, n);
}

double** jacobi(double** lnorm, int n) {
    int i, j;
    // init v_mat
    double** v_mat;
    v_mat = matrix_allocation(n, n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                v_mat[i][i] = 1.0;
            }
            else {
                v_mat[i][j] = 0.0;
            }
        }
    }
    double a_off = 0.0;
    for (i = 0; i < n; i++) {
        for (j = i+1; j < n; j++) {
            a_off += 2*pow(lnorm[i][j], 2.0);
        }
    }
    int current_iter = 0;
    while (current_iter < MAX_ROT) {
        p_rot_mat(lnorm, v_mat, &a_off, &current_iter, n);
        current_iter++;
    }
    // lnorm is changed in place
    return v_mat;
}


double** readfile(char* filename, int* n, int* d) {
    // read dimensions & number of points
    char c = 0;
    double num = 0;
    int point_i = 0;
    int i = 0;
    int dim = 0;
    int n_local = 0;
    int d_local = 0;
    int index = 0;
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("An Error Has Occurred\n");
        return 1;
    }
    while (fscanf(file, "%lf", &num) != EOF) {
        if ((c = getc(file)) == '\n') {
            (n_local)++;
            (d_local)++;
        }
        else if (c == ',') {
            (d_local)++;
        }
        else if (c == EOF) {
            (n_local)++;
            (d_local)++;
            break;
        }
        else {
            (n_local)++;
            (d_local)++;
            break;
        }
    }
    d_local /= n_local;
    fseek(file, 0L, SEEK_SET);
    // init points array
    double** allpoints = matrix_allocation(n_local, d_local);
    // read points from file
    while (fscanf(file, "%lf", &num) != EOF) {
        if ((c = getc(file)) == '\n') {
            allpoints[index][dim] = num;
            index++;
            dim = 0;
        }
        else if (c == ',') {
            allpoints[index][dim] = num;
            dim++;
        }
        else {
            allpoints[index][dim] = num;
            break;
        }
    }
    fclose(file);
    *n = n_local;
    *d = d_local;
    print_matrix(allpoints, *n, *d);
    return allpoints;
}


double** debug_a() {
    // init points array
    double** a_debug = matrix_allocation(3, 3);
    a_debug[0][0] = 3.0;
    a_debug[0][1] = 2.0;
    a_debug[0][2] = 4.0;
    a_debug[1][0] = 2.0;
    a_debug[1][1] = 0.0;
    a_debug[1][2] = 2.0;
    a_debug[2][0] = 4.0;
    a_debug[2][1] = 2.0;
    a_debug[2][2] = 3.0;
    return a_debug;
}
double** spkmeans(char* goal, char* filename, int* mat_size, int k) {
    int n =0, d = 0; int i, j;
    double** points = readfile(filename, mat_size, &d);
    n = *mat_size;
    double** w_mat = weighted_mat(points, n, d);
    if (strcmp(goal, "wam") == 0) {
        print_matrix(w_mat, n, n);
        free_matrix_allocation(points, n);
        return w_mat;
    }
    double** diag = diagonal_mat(n, w_mat);
    if (strcmp(goal, "ddg") == 0) {
        print_matrix(diag, n, n);
        free_matrix_allocation(points, n);
        free_matrix_allocation(w_mat, n);
        return diag;
    }
    double** lnorm = l_norm(n, w_mat,diag);
    if (strcmp(goal, "lnorm") == 0) {
        print_matrix(lnorm, n, n);
        free_matrix_allocation(points, n);
        free_matrix_allocation(w_mat, n);
        free_matrix_allocation(diag, n);
        return lnorm;
    }
    free_matrix_allocation(lnorm, n);
    double** v_mat;
    // L_Norm -> A (Eigenvalues of L_Norm)
    double** debuga = debug_a();
    v_mat = jacobi(debuga, 3);
    if (strcmp(goal, "jacobi") == 0) {
        print_matrix(lnorm, n, n);
        print_matrix(v_mat, n, n);
        free_matrix_allocation(points, n);
        free_matrix_allocation(w_mat, n);
        free_matrix_allocation(diag, n);
        return lnorm;
    }
}

int main(int argc, char* argv[]) {
    int n;
    if (argc != 3) {
        printf("Invalid input!");
        exit(1);
    }
    char* goal = argv[1];
    char* input_text = argv[2];
    if (strcmp(goal, "wam") == 0 || strcmp(goal, "ddg") == 0 || strcmp(goal, "lnorm") == 0 || strcmp(goal, "jacobi") == 0) {
        double** mat = spkmeans(goal, input_text, &n, 0);
        free_matrix_allocation(mat, n);
    }
}