#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

int main(int argc, char *argv[]){
    int d=0;
    int n=0;
    int max_iter=200;
    char check_string[100];
    char* input_text;
    char* output_text;
    char c;
    double num;
    int index = 0;
    int dim = 0;
    FILE* file;
    FILE* out;
    int less_than_eps = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int p = 0;
    int point_i = 0;
    double epsilon = 0.001 * 0.001;
    double** prev_centers;
    double** centers;
    double** buckets;
    double* euclidean_norm;
    double** allpoints;
    int* amounts;
    int K;

    if (argc == 4) {
        if (sscanf(argv[1],"%d",&K)!=1 || sscanf(argv[2], "%s", check_string)!=1 || sscanf(argv[3], "%s", check_string)!=1){
            printf("Invalid Input!\n");
            exit(1);
        }
        input_text = argv[2];
        output_text = argv[3];
    }
    else if (argc == 5) {
        if (sscanf(argv[1],"%d",&K)!=1 || sscanf(argv[2],"%d",&max_iter)!=1 || sscanf(argv[3], "%s", check_string)!=1 || sscanf(argv[4], "%s", check_string)!=1){
            printf("Invalid Input!\n");
            exit(1);
        }
        max_iter = atoi(argv[2]);
        input_text = argv[3];
        output_text = argv[4];
    }
    else {
        printf("Invalid Input!\n");
        exit(1);
    }
   
    file = fopen(input_text, "r");
    if (file == NULL){
        printf("An Error Has Occurred\n");
        return 1;
    }

    while (fscanf(file, "%lf", &num) != EOF){
        if ((c=getc(file)) == '\n'){
            n++;
            d++;
        }
        else if (c == ','){
            d++;
        }
        else if (c == EOF){
            n++;
            d++;
            break;
        }
        else{
            n++;
            d++;
            break;
        }
    }
    d /= n;
    
    fseek(file, 0L, SEEK_SET);
    amounts = calloc(K, sizeof(int));
    euclidean_norm = calloc(K, sizeof(double));    
    if (!amounts || !euclidean_norm) printf("problem with amounts/euclidean_norm");

    buckets = calloc(K, sizeof(double*));
    centers = calloc(K, sizeof(double*));
    prev_centers = calloc(K, sizeof(double*));
    if (!buckets || !centers ||!prev_centers) printf("problem with buckets/centers/prev_centers");
    for (p=0; p<K; p++){
        buckets[p] = calloc(d, sizeof(double));
        centers[p] = calloc(d, sizeof(double));
        prev_centers[p] = calloc(d, sizeof(double));
    }

    allpoints = calloc(n, sizeof(double*));
    if (!allpoints) printf("problem with allpoints");
    for (point_i=0; point_i <n; point_i++){
        allpoints[point_i] = calloc(d, sizeof(double));
    }

    while (fscanf(file, "%lf", &num) != EOF){
        if ((c=getc(file)) == '\n'){
            allpoints[index][dim] = num;
            index++;
            dim=0;   
        }
        else if (c == ','){
            allpoints[index][dim] = num;
            dim++;
        }
        else {
            allpoints[index][dim] = num;
            break;
        }
    }
    for (i=0; i<K; i++){
        for (j=0; j<d;j++){
            centers[i][j] = allpoints[i][j];
        }
    }
 
    for ( i = 0; i<max_iter; i++){
        for (p=0; p<n; p++){
            double min_dist = DBL_MAX;
            int cluster=0;
            double* point = allpoints[p];
            for(k = 0; k<K; k++){
                double dist = 0.0;
                for(j=0; j<d; j++){
                    dist += pow(centers[k][j]-point[j], 2.0);
                }
                if (dist < min_dist){
                    cluster = k;
                    min_dist = dist;
                }
            }
            amounts[cluster]++;
            for (j = 0; j<d; j++){
                buckets[cluster][j] += point[j];
            }
        }
        less_than_eps = 0;
        for (j=0; j<K; j++){
            for (k=0; k<d; k++){
                prev_centers[j][k] = centers[j][k];
                centers[j][k] = buckets[j][k]/(double)amounts[j];
                euclidean_norm[j] += pow(centers[j][k] - prev_centers[j][k], 2.0);
                buckets[j][k]=0.0;
            }
            amounts[j]=0;
            if (euclidean_norm[j] < epsilon) {
                less_than_eps += 1;
            }
        }

        if (less_than_eps == K) {
            break;
        }
    }
    fclose(file);
    
    out = fopen(output_text, "w");
    if (out == NULL){
        printf("An Error Has Occurred\n");
        return 1;
    }
    for (i=0; i<K*d; i++){
        if ( (i%d) != (d-1) ){
            fprintf(out, "%0.4f,", centers[(int)i/d][i%d]); 
        }
        else{
            fprintf(out, "%0.4f\n", centers[(int)i/d][i%d]);
        }
    }
    fclose(out);
 
    for (i = 0; i < n; i++) {
        free(allpoints[i]);
    }
    free(allpoints);
    for (p = 0; p < K; p++) {
        free(buckets[p]);
        free(centers[p]);
        free(prev_centers[p]);
    }
    free(buckets);
    free(centers);
    free(prev_centers);
   
    free(euclidean_norm);
    free(amounts);


    return 0;
}
