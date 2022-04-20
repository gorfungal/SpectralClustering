#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

static PyObject* spkmeans(PyObject *self, PyObject *args){
    char* filename = new char[100];
    int i;
    char* goal = new char[6];
    if(!PyArg_ParseTuple(args, "iss", &k, goal, filename)){
        return NULL;
    }
    double* mat = spkmean(goal, filename, k);
    int k = mat[0];
    int n = mat[1];
    PyObject* py_mat;
    py_mat = PyList_New((Py_ssize_t) (2+ (k*n)));
    PyList_SET_ITEM(py_mat, (Py_ssize_t) 0, PyFloat_FromDouble((double)k));
	for (i=0; i<k*n; i++){
		PyList_SET_ITEM(py_mat, (Py_ssize_t) (i+2), PyFloat_FromDouble(mat[i+2]));
	}
    return py_mat;
}

static PyMethodDef spkMethods[] = {
    {"spkmeans", (PyCFunction)spkmeans, METH_VARARGS, PyDoc_STR("One function to rule them all.")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef spkModule = {
    PyModuleDef_HEAD_INIT,
    "spkModule",
    "SPKMeans",
    -1,
    spkMethods
};

PyMODINIT_FUNC
PyInit_spkModule(void){
    PyObject *m;
    m = PyModule_Create(&spkModule);
    if(!m){
        return NULL;
    }
    return m;
}