// MPI_MatMultiVec.cpp: 主项目文件。
// 作者：谢家昱 西安交通大学-热流科学与工程教育部重点实验室

#include "stdafx.h"
#include <iostream>
#include <istream>
#include <fstream>
#include "mpi.h"
#define MAXN 100
#define MAXM 100
#define INFNAME "Input.dat"
#define OUFNAME "Output.dat"

int main(int argc, char *argv[]){
	int n, m, rank, nthread;
	double mat[MAXN][MAXM], vec[MAXM], sol[MAXN];
	void InputData(int *n, int *m, double mat[][MAXM], double vec[]);
	void MatMultiVec(int s, int e, int m, double mat[][MAXM], double vec[], double sol[]);
	void OutputData(int n, double sol[], int nthread, double dt);

	InputData(&n, &m, mat, vec);
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nthread);
	double dt = MPI_Wtime();
	int nrow = n / nthread;
	nrow += (n % nthread == 0) ? 0 : 1;
	int startrow = rank * nrow;
	int endrow = (rank + 1) * nrow - 1;
	if (endrow > n - 1) endrow = n - 1;
	MatMultiVec(startrow, endrow, m, mat, vec, sol);
	if (rank != 0){
		MPI_Send(&startrow, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&endrow, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD);
		for (int i = 0; i < endrow - startrow + 1; i++){
			MPI_Send(&sol[i], 1, MPI_DOUBLE, 0, i + 2, MPI_COMM_WORLD);
		}
	}
	else{
		for (int source = 1; source < nthread; source++){
			MPI_Recv(&startrow, 1, MPI_INTEGER, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&endrow, 1, MPI_INTEGER, source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int i = startrow; i <= endrow; i++){
				MPI_Recv(&sol[i], 1, MPI_DOUBLE, source, i - startrow + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		dt -= MPI_Wtime();
		OutputData(n, sol, nthread, -dt);
	}
	MPI_Finalize();
	return 0;
}

void InputData(int *n, int *m, double mat[][MAXM], double vec[]){
	using namespace std;
	ifstream fin(INFNAME);
	fin >> *n >> *m;
	for (int i = 0; i < *n; i++)
		for (int j = 0; j < *m; j++){
			fin >> mat[i][j];
		}
	for (int j = 0; j < *m; j++){
		fin >> vec[j];
	}
	fin.close();
	return;
}

void MatMultiVec(int s, int e, int m, double mat[][MAXM], double vec[], double sol[]){
	if (e - s + 1 <= 0) return;
	for (int i = s; i <= e; i++){
		double sum = 0.;
		for (int j = 0; j < m; j++){
			sum += mat[i][j] * vec[j];
		}
		sol[i - s] = sum;
	}
	return;
}

void OutputData(int n, double sol[], int nthread, double dt){
	using namespace std;
	ofstream fout(OUFNAME);
	fout << "Calculation Took " << dt << "s with " << nthread << " Threads.\n" << endl;
	fout << "Solution Vector:" << endl;
	for (int i = 0; i < n; fout << sol[i++] << endl);
	fout.close();
	return;
}