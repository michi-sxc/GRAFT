#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "sdust.h"

static PyObject* dust_mask(PyObject* self, PyObject* args) {
    const char* sequence;
    int length, T = 20, W = 64;
    
    if (!PyArg_ParseTuple(args, "s#|ii", &sequence, &length, &T, &W))
        return NULL;
    
    uint64_t* result = sdust(NULL, (const uint8_t*)sequence, length, T, W, NULL);
    
    if (result == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "DUST masking failed");
        return NULL;
    }
    
    PyObject* masked_seq = PyUnicode_FromStringAndSize(sequence, length);
    if (masked_seq == NULL) {
        free(result);
        return NULL;
    }
    
    Py_ssize_t i = 0;
    while (result[i] != (uint64_t)-1) {
        uint32_t start = (uint32_t)(result[i] >> 32);
        uint32_t end = (uint32_t)result[i];
        for (uint32_t j = start; j < end; j++) {
            PyUnicode_WriteChar(masked_seq, j, 'N');
        }
        i++;
    }
    
    free(result);
    return masked_seq;
}

static PyMethodDef DustMethods[] = {
    {"dust_mask", dust_mask, METH_VARARGS, "Apply DUST masking to a sequence."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef dustmodule = {
    PyModuleDef_HEAD_INIT,
    "dust_module",
    NULL,
    -1,
    DustMethods
};

PyMODINIT_FUNC PyInit_dust_module(void) {
    return PyModule_Create(&dustmodule);
}