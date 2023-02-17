#include <iostream>
#include "../header files/visual.h"
#include "Python.h"

void trisurf() {
	Py_Initialize();

	PyObject* obj = Py_BuildValue("s", "visualizers/trisurf.py");
	FILE* file = _Py_fopen_obj(obj, "r+");
	if (file != NULL) {
		PyRun_SimpleFile(file, "visualizers/trisurf.py");
	}
}

