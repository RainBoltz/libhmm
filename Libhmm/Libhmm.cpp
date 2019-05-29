// Libhmm.cpp : 定義 DLL 應用程式的匯出函式。
//

#include "stdafx.h"
#include <random>
#include <iostream>
#include "Libhmm.h"

namespace Libhmm {
	double* init_prob(int i) {
		double *M = new double[i];
		for (int c = 0; c < i; c++)
			M[c] = 1.0 / i;
		return M;
	}
	double** init_prob(int i, int j) {
		double **M = new double*[i];
		for (int b = 0; b < i; b++) {
			M[b] = new double[j];
			for (int c = 0; c < j; c++)
				M[b][c] = 1.0 / j;
		}
		return M;
	}
	void free_memory(double* M) {
		delete[] M;
		M = NULL;
	}
	void free_memory(double** M) {
		int nrow = _msize(M) / sizeof(**M);
		for (int i = 0; i < nrow; i++) {
			delete[] M[i];
		}
		delete[] M;
		M = NULL;
	}

	HMM::HMM(int n, int m): N(n), M(m){
		this->A = init_prob(n, n);
		this->B = init_prob(n, m);
		this->pi = init_prob(n);
	}
	HMM::HMM(char* model_path){
		load_model(model_path);
	}
	HMM::~HMM() {
		free_memory(this->A);
		free_memory(this->B);
		free_memory(this->pi);
	};
	void HMM::load_model(char* model_path){

	}
	void HMM::dump_model(char* output_path){

	}



	HMM_Trainer::HMM_Trainer(int n, int m, int t){

	}
	HMM_Trainer::~HMM_Trainer(){

	}

}