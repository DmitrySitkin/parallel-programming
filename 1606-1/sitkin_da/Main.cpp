//#include"mpi.h"
//#include<iostream>
//#include<time.h>
//using namespace std;
//int* CreateRandVec(int size)
//{
//	int *mas = new int[size];
//	for (int i = 0; i<size; i++)
//	{
//		mas[i] = rand();
//		//cout << vec[i] << endl;
//	}
//	return mas;
//}
////int FindMax(int* mas, int size, int _ProcNum, int _ProcRank)
////{
////	/*int _max=0;
////	for(int i=0;i<sendcounts[ProcRank];i++)
////	return _max;*/
////}
//int main(int argc, char *argv[])
//{
//	int n, ProcNum, ProcRank, *buf,_max;
//	int vecsize = atoi (argv[1]);
//	//cout << argv[0] << argv[1];// << argv[2] << argv[3];
//	double time = 0.0;
//	MPI_Status status;
//	srand(1);
//	//vec=CreateRandVec(vecsize);
//	//vec[9999] = 999999999;
//	MPI_Init(&argc, &argv);
//	//time = MPI_Wtime();
//	//cout << "2" << endl;
//	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
//	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
//	cout<<"proc" << ProcRank << endl;
//	MPI_Bcast(&vecsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
//	int *vec=nullptr;
//	//cout << "Bcast-" <<ProcRank<< endl;
//	if (ProcRank == 0)
//	{
//		vec = CreateRandVec(vecsize);
//		vec[9999] = 999999999;
//		
//	}
//	//else vec = nullptr;
//	cout << ProcRank << endl;
//	int b = vecsize%ProcNum;
//	int *sendcounts = new int[ProcNum]; 
//	int *displs = new int[ProcNum];
//	for (int i = 0; i < ProcNum; i++)
//	{
//		sendcounts[i] = vecsize / ProcNum;
//		if (b)
//		{
//			sendcounts[i]++;
//			b--;
//		}
//		displs[i] = sendcounts[i] * i;
//	}
//	int *recvbuf = new int[sendcounts[ProcRank]];
//	MPI_Scatterv(vec, sendcounts, displs, MPI_INT, recvbuf, sendcounts[ProcRank], MPI_INT, 0, MPI_COMM_WORLD);
//	/*if (ProcRank == 0)
//	{
//		for (int i = 1; i < ProcNum; i++)
//		{
//			MPI_Send(&vec[vecsize / ProcNum*i], vecsize / ProcNum, MPI_INT, i, 0, MPI_COMM_WORLD);
//			cout << "send-" << ProcRank << "->" << i << endl;
//		}
//	}
//	else
//	{
//		MPI_Recv(&vec[vecsize / (ProcNum*ProcRank)], vecsize / ProcNum, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
//		cout << "recv-" << ProcRank << "<-0" << endl;
//		
//	}*/
//	_max = vec[0];
//	for (int i = 0; i < sendcounts[ProcRank]; i++)
//	{
//		if (recvbuf[i] > _max) _max = recvbuf[i];
//	}
//	cout  <<"working-"<< ProcRank<<_max << endl;
//	/*if ((vecsize%ProcNum)&&(ProcRank==ProcNum-1))
//	{
//		for (int i = (vecsize / ProcNum)*ProcRank; i < vecsize; i++)
//		{
//			if (vec[i] > _max) _max = vec[i];
//		}
//	}
//	for (int i = (vecsize/ProcNum)*ProcRank; i < (vecsize / ProcNum)*(ProcRank+1); i++)
//	{
//		if (vec[i] > _max) _max = vec[i];
//	}*/
//	//_max = FindMax(vec, vecsize, ProcNum, ProcRank);
//	cout << ProcRank << endl;
//	int max = 0;
//	MPI_Barrier(MPI_COMM_WORLD);
//	MPI_Reduce(&_max, &max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
//	cout << ProcRank << endl;
////	time = MPI_Wtime() - time;
//	MPI_Finalize();
//	cout <<"max="<< max << endl;
//	//cout <<ProcRank << "_max=" << _max << endl;
//	cout << "time=" << time << endl;
//	return 0;
//}
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