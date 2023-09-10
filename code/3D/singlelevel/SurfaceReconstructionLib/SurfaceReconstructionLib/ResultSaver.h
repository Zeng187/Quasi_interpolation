#pragma once
#include <mat.h>
#include <string>

class ResultSaver {
public:
	ResultSaver();
	void Save(int res, double* result, const std::string& name);
	std::string description;
};

ResultSaver::ResultSaver() {
	description = "";
}

void ResultSaver::Save(int res, double* result, const std::string& name) {
	MATFile *pmatFile = NULL;
	mxArray *pMxArray = NULL;
	pmatFile = matOpen(("results/" + name + ".mat").c_str(), "w");
	int size = (res + 1)*(res + 1)*(res + 1);
	pMxArray = mxCreateDoubleMatrix(size, 1, mxREAL);
	memcpy((void*)mxGetPr(pMxArray), (void*)result, size * sizeof(double));
	matPutVariable(pmatFile, "V", pMxArray);
	mxDestroyArray(pMxArray);
	matClose(pmatFile);
}

