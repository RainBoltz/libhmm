﻿#pragma once

namespace Libhmm
{
	double* init_prob(int i);
	double** init_prob(int i, int j);
	void free_memory(double* M);
	void free_memory(double** M);

	public ref class HMM
	{
	public:
		HMM(int n, int m);
		HMM(char* model_path);
		~HMM();
		void load_model(char* model_path);
		void dump_model(char* output_path);

	private:
		int N;		/* number of hidden states;  Q={1,2,...,N} */
		int M;		/* number of observation symbols; V={1,2,...,M} */
		int T;		/* max path length */
		double** A;	/* A[1..N][1..N]. a[i][j] is the transition prob of going from state i at time t to state j at time t+1 */
		double** B;	/* B[1..N][1..M]. b[j][k] is the probability of of observing symbol k in state j */
		double* pi;	/* pi[1..N] pi[i] is the initial state distribution */
	};

	public ref class HMM_Trainer
	{
	public:
		HMM_Trainer(int n, int m, int t);
		~HMM_Trainer();
		//double **alpha[T][N], **beta[T][N];        /* evaluation problem */
		//double **delta[T][N]; int **phi[T][N];    /* decoding problem */
		//double **gamma[T][N], ***xi[T][N][N];     /* learning problem */
	};

}