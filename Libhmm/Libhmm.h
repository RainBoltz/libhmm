#pragma once

namespace Libhmm
{
	double* init_prob(int i);
	double** init_prob(int i, int j);
	double*** init_prob(int i, int j, int k);
	int** init_state(int i, int j);
	void free_memory(double* M);
	void free_memory(double** M);
	void free_memory(double*** M);

	public class HMM
	{
	public:
		HMM();
		HMM(int n, int m);
		HMM(char* model_path);
		~HMM();
		void init(int n, int m);
		void load_model(char* model_path);
		void dump_model(char* output_path);
		double decode_prob(int* o, double** delta, int** phi);	/* Viterbi Algorithm */
		int* decode_path(int* o, double** delta, int** phi);	/* Viterbi Algorithm */
		void learn(int* o, double** alpha, double** beta, double **gamma, double ***xi);	/* EM Algorithm */

	private:
		int N;		/* number of hidden states;  Q={1,2,...,N} */
		int M;		/* number of observation symbols; V={1,2,...,M} */
		int T;		/* max path length */
		double** A;	/* A[1..N][1..N]. a[i][j] is the transition prob of going from state i at time t to state j at time t+1 */
		double** B;	/* B[1..N][1..M]. b[j][k] is the probability of of observing symbol k in state j */
		double* pi;	/* pi[1..N] pi[i] is the initial state distribution */
		double forward(int* o, int T, double** alpha); 	/* evaluation problem solution-1: forward method */
		double backward(int* o, int T, double** beta);	/* evaluation problem solution-2: backward method */
		
	};

	public ref class HMM_Trainer
	{
	public:
		HMM_Trainer(int n, int m, int t);
		~HMM_Trainer();
		double **alpha;	/* evaluation problem: forward (T >= 2)*/
		double **beta;	/* evaluation problem: backward (T >= 2)*/
		double **delta;	/* decoding problem: pi*B (T >= 2)*/
		int **phi;    	/* decoding problem: states */
		double **gamma;	/* learning problem variables */
		double ***xi;	/* learning problem variables */
		void train(int*** training_set, int nepoch);
		int recognize(int* o);
	private:
		bool init;
		int N, M, MaxT;
		HMM* hmms;
	};

}