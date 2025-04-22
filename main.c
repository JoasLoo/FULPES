#include <Python.h>

int main() {
    Py_Initialize();

    // Load script.py
    PyObject *pName = PyUnicode_DecodeFSDefault("script");
    PyObject *pModule = PyImport_Import(pName);
    Py_DECREF(pName);

    if (pModule != NULL) {
        // Call say_hello()
        PyObject *pFunc = PyObject_GetAttrString(pModule, "say_hello");
        if (PyCallable_Check(pFunc)) {
            PyObject_CallObject(pFunc, NULL);
        }
        Py_XDECREF(pFunc);

        // Call add(3, 5)
        pFunc = PyObject_GetAttrString(pModule, "add");
        if (PyCallable_Check(pFunc)) {
            PyObject *args = PyTuple_Pack(2, PyLong_FromLong(3), PyLong_FromLong(5));
            PyObject *pValue = PyObject_CallObject(pFunc, args);
            Py_DECREF(args);
            if (pValue != NULL) {
                long result = PyLong_AsLong(pValue);
                printf("Result of add(3, 5): %ld\n", result);
                Py_DECREF(pValue);
            }
        }
        Py_XDECREF(pFunc);

        Py_DECREF(pModule);
    } else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"script.py\"\n");
    }

    Py_Finalize();
    return 0;
}
