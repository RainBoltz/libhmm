// Libhmm.cpp : 定義 DLL 應用程式的匯出函式。
//

#include "stdafx.h"
#include <random>
#include <iostream>
#include "Libhmm.h"

namespace Libhmm {
	void init_prob(double* M, int i) {
		M = new double[i];
		for (int c = 0; c < i; c++)
			M[c] = 1.0 / i;
	}
	void init_prob(double** M, int i, int j) {
		for (int a = 0; a < i; i++) {
			M = new double*[i];
			for (int b = 0; b < j; b++) {
				M[b] = new double[j];
				for (int c = 0; c < j; c++)
					M[b][c] = 1.0 / j;
			}
		}
	}
	void free_memory(double* M) {
		delete[] M;
	}
	void free_memory(double** M) {
		int _n = _msize(M[0]) / sizeof(M[0][0]);
		for (int i = _n - 1; i >= 0; i++)
			delete[] M[i];
	}

	HMM::HMM(int n, int m): N(n), M(m){
		init_prob(A, n, n);
		init_prob(B, n, m);
		init_prob(pi, n);
		std::cout << "CREATED!" << std::endl;
	}
	HMM::HMM(char* model_path){
		load_model(model_path);
	}
	HMM::~HMM() {
		free_memory(A);
		free_memory(B);
		free_memory(pi);
		std::cout << "DELETED!" << std::endl;
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