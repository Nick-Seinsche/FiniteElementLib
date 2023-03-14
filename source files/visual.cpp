#include <iostream>
#include "Python.h"
#include <cstdlib>

#include "../header files/visual.h"


const wchar_t* GetWC(const char* c){
	const size_t cSize = strlen(c) + 1;
	size_t size = strlen(c) + 1;
	wchar_t* portName = new wchar_t[size];
	size_t outSize;
	mbstowcs_s(&outSize, portName, size, c, size-1);
	return portName;
}

wchar_t** GetWCStar(char** c, int sz) {
	wchar_t** av = new wchar_t* [sz];
	for (int i = 0; i < sz; i++) {
		av[i] = (wchar_t*) GetWC(c[i]);
	}
	return av;
}

void visual::trisurf(std::string name_triang, std::string name_sol) {
	int argc = 5;
	char* argv[5];

	argv[0] = (char*) "trisurf.py";
	argv[1] = (char*) "-t";
	argv[2] = (char*) name_triang.c_str();
	argv[3] = (char*) "-s";
	argv[4] = (char*) name_sol.c_str();


	Py_SetProgramName(GetWC(argv[0]));
	Py_Initialize();
	PySys_SetArgv(argc, GetWCStar(argv, argc));

	PyObject* obj = Py_BuildValue("s", "visualizers/trisurf.py");
	FILE* file = _Py_fopen_obj(obj, "r+");
	if (file != NULL) {
		PyRun_SimpleFile(file, "visualizers/trisurf.py");
		Py_Finalize();
	}
}

void visual::tricontour(std::string name_triang, std::string name_sol) {
	int argc = 5;
	char* argv[5];

	argv[0] = (char*)"tricontour.py";
	argv[1] = (char*)"-t";
	argv[2] = (char*)name_triang.c_str();
	argv[3] = (char*)"-s";
	argv[4] = (char*)name_sol.c_str();


	Py_SetProgramName(GetWC(argv[0]));
	Py_Initialize();
	PySys_SetArgv(argc, GetWCStar(argv, argc));

	PyObject* obj = Py_BuildValue("s", "visualizers/tricontour.py");
	FILE* file = _Py_fopen_obj(obj, "r+");
	if (file != NULL) {
		PyRun_SimpleFile(file, "visualizers/tricontour.py");
		Py_Finalize();
	}
}
