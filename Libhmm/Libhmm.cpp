// Libhmm.cpp : 定義 DLL 應用程式的匯出函式。
//

#include "stdafx.h"
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
	}
	void HMM::load_model(char* model_path){
		//TODO
	}
	void HMM::dump_model(char* output_path){
		//TODO
	}

	double HMM::forward(int* o, int T, double** alpha){
		//int T = _msize(o) / sizeof(*o);
		for (int t=0; t<T; ++t)
			for (int j=0; j<this->N; ++j)
				if (t == 0)
					alpha[t][j] = this->pi[j] * this->B[j][o[t]];
				else
				{
					double p = 0;
					for (int i=0; i<this->N; ++i)
						p += alpha[t-1][i] * this->A[i][j];
					alpha[t][j] = p * this->B[j][o[t]];
				}
		double p = 0;
		for (int i=0; i<this->N; ++i)
			p += alpha[T-1][i];
		return p;
	}
	double HMM::backward(int* o, int T, double** beta){
		//int T = _msize(o) / sizeof(*o);
		for (int t=T-1; t>=0; --t)
			for (int i=0; i<this->N; ++i)
				if (t == T-1)
					beta[t][i] = 1.0;
				else
				{
					double p = 0;
					for (int j=0; j<this->N; ++j)
						p += this->A[i][j] * this->B[j][o[t+1]] * beta[t+1][j];
					beta[t][i] = p;
				}
	
		double p = 0;
		for (int j=0; j<this->N; ++j)
			p += this->pi[j] * this->B[j][o[0]] * beta[0][j];
		return p;
	}
	double HMM::decode_prob(int* o, double** delta, double** phi){
		int T = _msize(o) / sizeof(*o);
		for (int t=0; t<T; ++t)
			for (int j=0; j<this->N; ++j)
				if (t == 0)
					delta[t][j] = this->pi[j] * this->B[j][o[t]];
				else
				{
					double p = -1e9;
					for (int i=0; i<this->N; ++i)
					{
						double w = delta[t-1][i] * this->A[i][j];
						if (w > p) p = w, phi[t][j] = i;
					}
					delta[t][j] = p * this->B[j][o[t]];
				}
	
		double p = -1e9;
		for (int j=0; j<this->N; ++j)
			if (delta[T-1][j] > p)
				p = delta[T-1][j], q[T-1] = j;
		
		return p;
	}
	int* HMM::decode_path(int* o, double** delta, double** phi){
		//q: best sequence of hidden states for observation o
		int T = _msize(o) / sizeof(*o);
		int *q = new int[T];
		for (int t=0; t<T; ++t)
			for (int j=0; j<this->N; ++j)
				if (t == 0)
					delta[t][j] = this->pi[j] * this->B[j][o[t]];
				else
				{
					double p = -1e9;
					for (int i=0; i<this->N; ++i)
					{
						double w = delta[t-1][i] * this->A[i][j];
						if (w > p) p = w, phi[t][j] = i;
					}
					delta[t][j] = p * this->B[j][o[t]];
				}
		
		for (int t=T-1; t>0; --t)
			q[t-1] = phi[t][q[t]];
		
		return q;
	}
	void HMM::learn(int* o, double** alpha, double** beta, double **gamma, double ***xi){
		int T = _msize(o) / sizeof(*o);
		this->forward(o, T);
		this->backward(o, T);
	
		for (int t=0; t<T; ++t)
		{
			double p = 0;
			for (int i=0; i<this->N; ++i)
				p += alpha[t][i] * beta[t][i];
			assert(p != 0);
	
			for (int i=0; i<this->N; ++i)
				gamma[t][i] = alpha[t][i] * beta[t][i] / p;
		}
	
		for (int t=0; t<T-1; ++t)
		{
			double p = 0;
			for (int i=0; i<this->N; ++i)
				for (int j=0; j<this->N; ++j)
					p += alpha[t][i] * this->A[i][j] * this->B[j][o[t+1]] * beta[t+1][j];
			assert(p != 0);
	
			for (int i=0; i<this->N; ++i)
				for (int j=0; j<this->N; ++j)
					xi[t][i][j] = alpha[t][i] * this->A[i][j] * this->B[j][o[t+1]] * beta[t+1][j] / p;
		}
	
		// update pi
		for (int i=0; i<this->N; ++i)
			 this->pi[i] = gamma[0][i];

		// update A
		for (int i=0; i<this->N; ++i)
		{
			double p2 = 0;
			for (int t=0; t<T-1; ++t)
				p2 += gamma[t][i];
			assert(p2 != 0);

			for (int j=0; j<this->N; ++j)
			{
				double p1 = 0;
				for (int t=0; t<T-1; ++t)
					p1 += xi[t][i][j];
				 this->A[i][j] = p1 / p2;
			}
		}
	
		// update B
		for (int i=0; i<this->N; ++i)
		{
			double p[this->M] = {0}, p2 = 0;
			for (int t=0; t<T; ++t)
			{
				p[o[t]] += gamma[t][i];
				p2 += gamma[t][i];
			}
			assert(p2 != 0);
	
			for (int k=0; k<this->M; ++k)
				this->B[i][k] = p[k] / p2;
		}
	}

	HMM_Trainer::HMM_Trainer(int n, int t){
		// t is the max limit of observation length
		this->alpha = init_prob(t, n);
		this->beta = init_prob(t, n);
		this->delta = init_prob(t, n);
		this->phi = init_state(t, n);
		this->gamma = init_prob(t, n);
		this->xi = init_prob(t, n, n);
	}
	HMM_Trainer::~HMM_Trainer(){
		freeMemory(this->alpha);
		freeMemory(this->beta);
		freeMemory(this->delta);
		freeMemory(this->phi);
		freeMemory(this->gamma);
		freeMemory(this->xi);
	}
	void HMM_Trainer::train(HMM::HMM* hmms, int*** training_set, int nepoch){
		int nclass = _msize(training_set) / sizeof(***training_set);
		int ndata = _msize(*training_set) / sizeof(***training_set);
		for (int i=0; i<nclass; ++i)
			for (int l=0; l<nepoch; ++l)
				for (int j=0; j<ndata; ++j){
					if(l==0) std::cout << "START TRAINING" << std::endl;
					else if(l%10==0) std::cout << "TRAINED " << l+1 << "EPOCHS..." << std::endl;
					hmms[i]->learn(training_set[i], this->alpha, this->beta, this->gamma, this->xi);
				}
		std::cout << "TRAINING FINISHED!" << std::endl;
	}
	int classify(HMM::HMM* hmms, int* o){
		int T = _msize(o) / sizeof(*o);
		int nclass = _msize(hmms) / sizeof(*hmms);
		int ans = -1;
		double p = -1e9;
		for (int i=0; i<nclass; ++i)
		{
			double pp = hmms[i]->decode_prob(o, T);
			if (pp > p) p = pp, ans = i;
		}
		return ans;
	}
	
}