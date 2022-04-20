#ifndef SPKMEANS_H
#define SPKMEANS_H

double** mat_mul(int n, double** m1, double** m2);
double** weighted_mat(double** allpoints, int n, int d);
double** diagonal_mat(int n, double** weighted);
double** l_norm(int n, double** weighted);
void p_rot_mat(double** a_mat,double** v_mat,double* a_off_diag_sum, int* current_iter, int n);
int compare(float num1, float num2);
int num_of_clusters(float **eigenvalues);
double** jacobi(double** lnorm, int n);
int Eigengap_heuristic(double** a_mat, double** v_mat, int n);
double** spkmean(char* goal, char* filename, int k);

#endif