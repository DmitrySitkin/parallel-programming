
#include <iostream>
#include "mpi.h"
#include <ctime>
#include <limits>
using namespace std;
void FillingArray(double array[], int n) 
{
	srand(time(nullptr));
	for (int i = 0; i < n; i++)
	{
		double ri = (double)rand() *(double)rand()/10000;
		array[i] = ri;
	}
}
int main(int argc, char **argv)
{
	int vecsize = 10;
	int ProcNum, ProcRank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) 
	{
		vecsize = atoi(argv[1]);
	}

	MPI_Bcast(&vecsize, 1, MPI_INT, 0, MPI_COMM_WORLD); 
	double tstart, tfinish;
	double *vec;
	double totalMax = 0.0; 

	if (ProcRank == 0)
	{
		vec = new double[vecsize];
		FillingArray(vec, vecsize);
		tstart = MPI_Wtime();
	}
	else {
		vec = nullptr;
	}

	int a = vecsize % ProcNum;

	int *sendcounts = new int[ProcNum]; 
	int *displs = new int[ProcNum]; 
	int sum = 0; 
	for (int i = 0; i < ProcNum; i++)
	{
		sendcounts[i] = vecsize / ProcNum;
		if (a)
		{
			sendcounts[i]++;
			a--;
		}
		displs[i] = sum;
		sum += sendcounts[i];
	}

	double *recbuf = new double[sendcounts[ProcRank]]; 

	double max = 0.0; 

	MPI_Scatterv(vec, sendcounts, displs, MPI_DOUBLE, recbuf, sendcounts[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD); 

	for (int i = 0; i < sendcounts[ProcRank]; i++) 
	{
		if (recbuf[i] > max) max = recbuf[i];
	}

	MPI_Reduce(&max, &totalMax, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); //reduces partial amount on all processes to a total amount
	tfinish = MPI_Wtime();
	double tstart1;
	if (ProcRank == 0)
	{
		cout << "MAX=" << totalMax << endl;
	}
	MPI_Finalize();
	delete[] recbuf, sendcounts, displs, vec;
	return 0;
}