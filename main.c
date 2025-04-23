#include <Python.h>

int main() {
    Py_Initialize(); //Initialize the Python interpreter.

    // Load script.py
    PyObject *pName = PyUnicode_DecodeFSDefault("script");
    PyObject *pModule = PyImport_Import(pName);     //Import "script".py
    Py_DECREF(pName);                               //Free memory used by pName.

    if (pModule != NULL) {
        // Call say_hello()
        PyObject *pFunc_hello = PyObject_GetAttrString(pModule, "say_hello"); //get hello func
        PyObject *pFunc_add = PyObject_GetAttrString(pModule, "add"); //get add function from script.py
        PyObject *args = PyTuple_Pack(3, PyLong_FromLong(3), PyLong_FromLong(5), PyLong_FromLong(100));

        //if (PyCallable_Check(pFunc_hello)) {                                  //check if the func is gotten properly
            PyObject_CallObject(pFunc_hello, NULL);                           //call function with no input args
        //}
        
        if (PyCallable_Check(pFunc_add)) {
            PyObject *pValue = PyObject_CallObject(pFunc_add, args);
            //
            if (pValue != NULL) {
                long result = PyLong_AsLong(pValue);
                printf("Result of add(3, 5): %ld\n", result);
                PyTuple_SetItem(args, 2, PyLong_FromLong(result)); //start counting from 0, so we change the third variable
                Py_DECREF(pValue);
            }
        }
        
        for (int i = 0; i < 0; i++) {
            long result1 = PyLong_AsLong(PyObject_CallObject(pFunc_add, args));
            printf("Result of add(3, 5, 100): %ld\n", result1);
            PyTuple_SetItem(args, 2, PyLong_FromLong(result1));
        }

        //PyObject_CallObject(pFunc_hello, NULL); 
        //free up space
        Py_XDECREF(pFunc_hello);
        Py_XDECREF(pFunc_add);
        Py_XDECREF(args);
        Py_XDECREF(pModule);
    } else {
        PyErr_Print();
        fprintf(stderr, "Failed to load \"script.py\"\n");
    }

    Py_Finalize();  //exit the python environment
    return 0;
}
