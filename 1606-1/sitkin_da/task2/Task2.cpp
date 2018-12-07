#include"mpi.h"
#include<iostream>
#include <ctime>
using namespace std;
void my_bcast(void* data, 
	int count, 
	MPI_Datatype datatype, 
	int root,
	MPI_Comm communicator)
{
	int world_rank;
	MPI_Comm_rank(communicator, &world_rank);
	int world_size;
	MPI_Comm_size(communicator, &world_size);

	if (world_rank == root) {
		int i,k,m;
		int n = (int)log2(world_size);
		double mod = log2f(world_size);
		for(i=0;i<n;i++)
			for (k = 0; k < pow(2,i); k++)
			{
				//if (k != world_rank)
				
					//std::cout << k + pow(2, i)<<std::endl;
					MPI_Send(data, count, datatype, k+pow(2,i), 0, communicator);
				
			}
		if (world_size - pow(2, n) !=0)
		{
			m = world_size - pow(2, n);
				for (k =0; k < m; k++)
				{
					//if (k > n)

				//	std::cout << k + pow(2, n) << std::endl;
					MPI_Send(data, count, datatype, k+pow(2,n), 0, communicator);

				}
		}
	}
	else {
		MPI_Recv(data, count, datatype, root, 0, communicator,
			MPI_STATUS_IGNORE);
	}
}
void FillingArray(double array[], int n)
{
	srand(time(nullptr));
	for (int i = 0; i < n; i++)
	{
		double ri = (double)rand() *(double)rand() / 10000;
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
	double time, time1;
	if (ProcRank == 0)
	{
		vecsize = atoi(argv[1]);
	}
	if (ProcRank==0)
	{
		time = MPI_Wtime();
	}
	MPI_Bcast(&vecsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		time = MPI_Wtime()-time;
	}
	if (ProcRank == 0)
	{
		time1 = MPI_Wtime();
	}
	//my_bcast(&vecsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (ProcRank == 0)
	{
		time1 = MPI_Wtime() - time1;
	}
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

	MPI_Reduce(&max, &totalMax, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
	tfinish = MPI_Wtime();
	double tstart1;
	if (ProcRank == 0)
	{
		cout << "MPI_Bcast time=" << time << endl << "my_bcast time=" << time1 << endl;
		cout << "MAX=" << totalMax << endl;
	}
	MPI_Finalize();
	return 0;
}