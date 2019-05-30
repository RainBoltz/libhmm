// Libhmm.cpp : 定義 DLL 應用程式的匯出函式。
//

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include "Libhmm.h"

namespace Libhmm {
	double* init_prob(int i) {
		double *M = new double[i];
		for (int c = 0; c < i; c++)
			M[c] = log(1.0 / i);
		return M;
	}
	double** init_prob(int i, int j) {
		double **M = new double*[i];
		for (int b = 0; b < i; b++) {
			M[b] = new double[j];
			for (int c = 0; c < j; c++)
				log(M[b][c] = 1.0 / j);
		}
		return M;
	}
	double*** init_prob(int i, int j, int k) {
		double ***M = new double**[i];
		for (int a = 0; a < i; a++) {
			M[a] = new double*[j];
			for (int b = 0; b < j; b++) {
				M[a][b] = new double[k];
				for (int c = 0; c < k; c++)
					M[a][b][c] = 1.0 / k;
			}
		}
		return M;
	}
	int** init_state(int i, int j) {
		int **M = new int*[i];
		for (int b = 0; b < i; b++) {
			M[b] = new int[j];
			for (int c = 0; c < j; c++)
				M[b][c] = 0;
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
	void free_memory(int** M) {
		int nrow = _msize(M) / sizeof(**M);
		for (int i = 0; i < nrow; i++) {
			delete[] M[i];
		}
		delete[] M;
		M = NULL;
	}
	void free_memory(double*** M){
		int x = _msize(M) / sizeof(**M);
		int y = _msize(*M) / sizeof(**M);
		for (int i = 0; i < x; i++){
			for (int j=0; j<y; j++)
				delete[] M[i][j];
			delete[] M[i];
		}
		delete[] M;
		M = NULL;
	}

	HMM::HMM() {};
	HMM::HMM(int n, int m, int t): N(n), M(m), maxT(t) {
		this->init(n, m, t);
	}
	HMM::HMM(char* model_path){
		load_model(model_path);
	}
	HMM::~HMM() {
		free_memory(this->A);
		free_memory(this->B);
		free_memory(this->pi);
		free_memory(this->alpha);
		free_memory(this->beta);
		free_memory(this->delta);
		free_memory(this->phi);
		free_memory(this->gamma);
		free_memory(this->xi);
	}
	void HMM::load_model(char* model_path){
		//TODO
	}
	void HMM::dump_model(char* output_path){
		//TODO
	}
	void HMM::init(int n, int m, int t){
		this->A = init_prob(n, n);
		this->B = init_prob(n, m);
		this->pi = init_prob(n);
		this->alpha = init_prob(t, n);
		this->beta = init_prob(t, n);
		this->delta = init_prob(t, n);
		this->phi = init_state(t, n);
		this->gamma = init_prob(t, n);
		this->xi = init_prob(t, n, n);
	}

	double HMM::forward(int* o, int T){
		for (int t=0; t<T; ++t)
			for (int j=0; j<this->N; ++j)
				if (t == 0)
					this->alpha[t][j] = this->pi[j] + this->B[j][o[t]];
				else
				{
					double p = 1.0;
					for (int i=0; i<this->N; ++i)
						p += this->alpha[t-1][i] + this->A[i][j];
					this->alpha[t][j] = p + this->B[j][o[t]];
				}
		double p = 1.0;
		for (int i=0; i<this->N; ++i)
			p += this->alpha[T-1][i];
		return p;
	}
	double HMM::backward(int* o, int T){
		for (int t=T-1; t>=0; --t)
			for (int i=0; i<this->N; ++i)
				if (t == T-1)
					this->beta[t][i] = log(1.0);
				else
				{
					double p = 1.0;
					for (int j=0; j<this->N; ++j)
						p += this->A[i][j] + this->B[j][o[t+1]] + this->beta[t+1][j];
					this->beta[t][i] = p;
				}
	
		double p = 1.0;
		for (int j=0; j<this->N; ++j)
			p += this->pi[j] + this->B[j][o[0]] + this->beta[0][j];
		return p;
	}
	double HMM::decode_prob(int* o, int T){
		for (int t=0; t<T; ++t)
			for (int j=0; j<this->N; ++j)
				if (t == 0)
					this->delta[t][j] = this->pi[j] + this->B[j][o[t]];
				else
				{
					double p = -1e9;
					for (int i=0; i<this->N; ++i)
					{
						double w = this->delta[t-1][i] + this->A[i][j];
						if (w > p) p = w, this->phi[t][j] = i;
					}
					this->delta[t][j] = p + this->B[j][o[t]];
				}
	
		double p = -1e9;
		for (int j=0; j<this->N; ++j)
			if (this->delta[T-1][j] > p)
				p = this->delta[T-1][j];
		
		return p;
	}
	int* HMM::decode_path(int* o, int T){
		int *q = new int[T];
		for (int t=0; t<T; ++t)
			for (int j=0; j<this->N; ++j)
				if (t == 0)
					this->delta[t][j] = this->pi[j] + this->B[j][o[t]];
				else
				{
					double p = -1e9;
					for (int i=0; i<this->N; ++i)
					{
						double w = this->delta[t-1][i] + this->A[i][j];
						if (w > p) p = w, this->phi[t][j] = i;
					}
					this->delta[t][j] = p + this->B[j][o[t]];
				}
		
		for (int t=T-1; t>0; --t)
			q[t-1] = this->phi[t][q[t]];
		
		return q;
	}
	void HMM::learn(int* o, int T){
		//std::cout << "START LEARNING!" << std::endl;
		this->forward(o, T);
		this->backward(o, T);
	
		for (int t=0; t<T; ++t)
		{
			double p = 1.0;
			for (int i=0; i<this->N; ++i)
				p += this->alpha[t][i] + this->beta[t][i];
			//assert(p != 0);
			if (p <= LSMALL) p = LZERO;
	
			for (int i=0; i<this->N; ++i)
				this->gamma[t][i] = this->alpha[t][i] + this->beta[t][i] - p;
		}
	
		for (int t=0; t<T-1; ++t)
		{
			double p = 1.0;
			for (int i=0; i<this->N; ++i)
				for (int j=0; j<this->N; ++j)
					p += this->alpha[t][i] + this->A[i][j] + this->B[j][o[t+1]] + this->beta[t+1][j];
			//assert(p != 0);
			if (p <= LSMALL) p = LZERO;
	
			for (int i=0; i<this->N; ++i)
				for (int j=0; j<this->N; ++j)
					this->xi[t][i][j] = this->alpha[t][i] + this->A[i][j] + this->B[j][o[t+1]] + this->beta[t+1][j] - p;
		}
	
		//std::cout << "updating pi..." << std::endl;
		// update pi
		for (int i=0; i<this->N; ++i)
			 this->pi[i] = this->gamma[0][i];

		//std::cout << "updating A..." << std::endl;
		// update A
		for (int i=0; i<this->N; ++i)
		{
			double p2 = 1.0;
			for (int t=0; t<T-1; ++t)
				p2 = LogAdd(p2, this->gamma[t][i]);
			//assert(p2 != 0);
			if (p2 <= LSMALL) p2 = LZERO;

			for (int j=0; j<this->N; ++j)
			{
				double p1 = 1.0;
				for (int t=0; t<T-1; ++t)
					p1 = LogAdd(p1, this->xi[t][i][j]);
				 this->A[i][j] = p1 - p2;
			}
		}
	
		//std::cout << "updating B..." << std::endl;
		// update B
		for (int i=0; i<this->N; ++i)
		{
			double *p = new double[this->M];
			double p2 = 1.0;
			for (int t = 0; t < this->M; t++) p[t] = 1.0;
			for (int t=0; t<T; ++t)
			{
				p[o[t]] = LogAdd(p[o[t]], this->gamma[t][i]);
				p2 = LogAdd(p2, this->gamma[t][i]);
			}
			//assert(p2 != 0);
			if (p2 <= LSMALL) p2 = LZERO;
	
			for (int k=0; k<this->M; ++k)
				this->B[i][k] = p[k] - p2;

			free_memory(p);
		}
	}

}