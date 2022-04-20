#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


static PyObject* fit(int N, int lg, int D, int MAX, double e,double* Centers_init,double* Points){
    int d=D;
    int n=N;
    int K=lg;
    int max_iter=MAX;
    int less_than_eps = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int p = 0;
    int point_i = 0;
    double epsilon = e * e;
    double** prev_centers;
    double** centers;
	PyObject *centers_fin_object;
    double** buckets;
    double* euclidean_norm;
    double** allpoints;
    int* amounts;


    allpoints = calloc(n, sizeof(double*));
    if (!allpoints) printf("problem with memory allocation allpoints");
    for (point_i=0; point_i <n; point_i++){
        allpoints[point_i] = calloc(d, sizeof(double));
    }

    for (p=0; p<n; p++){
        for (i=0; i<d; i++){
            allpoints[p][i] = Points[p*d + i];
        }
    }

    amounts = calloc(K, sizeof(int));
    euclidean_norm = calloc(K, sizeof(double));
    if (!amounts || !euclidean_norm) printf("problem with memory allocation amounts/euclidean_norm\n");

    buckets = calloc(K, sizeof(double*));
    centers = calloc(K, sizeof(double*));

    prev_centers = calloc(K, sizeof(double*));
    if (!buckets || !centers ||!prev_centers) printf("problem with memory allocation buckets/centers/prev_centers");
    for (p=0; p<K; p++){
        buckets[p] = calloc(d, sizeof(double));
        centers[p] = calloc(d, sizeof(double));
        prev_centers[p] = calloc(d, sizeof(double));
    }

    for (p=0; p<K; p++){
        for (i=0; i<d; i++){
            centers[p][i] = Centers_init[p*d + i];
        }
    }

//start iterations
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
	
	centers_fin_object = PyList_New((Py_ssize_t) K*d);
	for (i=0; i<K; i++){
		for (j=0; j<d; j++){
			PyList_SET_ITEM(centers_fin_object,(Py_ssize_t) (i*d + j) , PyFloat_FromDouble(centers[i][j]));
		}
	}


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
    //free(centers);
    free(prev_centers);

    free(euclidean_norm);
    free(amounts);

    return centers_fin_object;
}


static PyObject* fit_i(PyObject *self, PyObject *args){
    int n;
    int k;
    int dim;
    int index = 0; 
    int max_iter;
    double epsilon;
    PyObject *centers_object;
    PyObject *points_object;
    double *points;
    double *centers_init;
    Py_ssize_t points_length;
    Py_ssize_t centers_length;
    PyObject *res;
    if(!PyArg_ParseTuple(args, "iiiidOO",&n, &k, &dim, &max_iter, &epsilon, &centers_object, &points_object )){
        return NULL;
    }
    //points to double
    points_length = PyObject_Length(points_object);
    if (points_length < 0)
        return NULL;
    points = (double *) malloc(sizeof(double *) * points_length);
    if (points == NULL)
        return NULL;
    for (index = 0; index < points_length; index++) {
        PyObject *item;
        item = PyList_GetItem(points_object, index);
        if (!PyFloat_Check(item)){
            return NULL;
        }
        points[index] = PyFloat_AsDouble(item);
    }
    //centers_init
    centers_length = PyObject_Length(centers_object);
    if (centers_length < 0)
        return NULL;
    centers_init = (double *) malloc(sizeof(double *) * centers_length);
    if (centers_init == NULL)
        return NULL;
    for (index = 0; index < centers_length; index++) {
        PyObject *item;
        item = PyList_GetItem(centers_object, index);
        if (!PyFloat_Check(item)){
            return NULL;
        }
        centers_init[index] = PyFloat_AsDouble(item);
    }
    res = Py_BuildValue("O", fit(n, k, dim, max_iter,epsilon,centers_init,points ));
    free(centers_init);
    free(points);

    return res;
    /* how to parse args with/out max_iter*/
}

static PyMethodDef kMethods[] = {
    {"fit", (PyCFunction) fit_i, METH_VARARGS, PyDoc_STR("info")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    kMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void){
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if(!m){
        return NULL;
    }
    return m;
}
