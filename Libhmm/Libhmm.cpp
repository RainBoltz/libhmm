// Libhmm.cpp : 定義 DLL 應用程式的匯出函式。
//

#include "stdafx.h"
#include <random>
#include "Libhmm.h"

namespace Libhmm {
	HMM::HMM(int n, int m): N(n), M(m){
		init_prob(A, n, n);
		init_prob(B, n, m);
		init_prob(pi, n);
	}
	HMM::HMM(char* model_path){
		load_model(model_path);
	}
	HMM::~HMM() {
		free_memory(A);
		free_memory(B);
		free_memory(pi);
	};
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
	void HMM::load_model(char* model_path){

	}
	void HMM::dump_model(char* output_path){

	}
}