#include"mpi.h"
#include<iostream>
#include<time.h>
using namespace std;
int* CreateRandVec(int size)
{
	int *mas = new int[size];
	for (int i = 0; i<size; i++)
	{
		mas[i] = rand();
		//cout << vec[i] << endl;
	}
	return mas;
}
int FindMax(int* mas, int size, int _ProcNum, int _ProcRank)
{
	int _max=0;
	if ((size%_ProcNum) && (_ProcRank == _ProcNum - 1))
	{
		for (int i = (size / _ProcNum)*_ProcRank; i < size; i++)
		{
			if (mas[i] > _max) _max = mas[i];
		}
	}
	for (int i = (size / _ProcNum)*_ProcRank; i < (size / _ProcNum)*(_ProcRank + 1); i++)
	{
		if (mas[i] > _max) _max = mas[i];
	}
	return _max;
}
int main(int argc, char *argv[])
{
	int n, ProcNum, ProcRank, *vec, *buf, max,_max;
	int vecsize = atoi (argv[1]);
	//cout << argv[0] << argv[1];// << argv[2] << argv[3];
	double time = 0.0;
	MPI_Status status;
	srand(1);
	vec=CreateRandVec(vecsize);
	vec[9999] = 999999999;
	MPI_Init(&argc, &argv);
	time = MPI_Wtime();
	//cout << "2" << endl;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	cout<<"proc" << ProcRank << endl;
	MPI_Bcast(&vecsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	cout << "Bcast-" <<ProcRank<< endl;
	if (ProcRank == 0)
	{
		vec = CreateRandVec(vecsize);
		vec[9999] = 999999999;
		for (int i = 1; i < ProcNum; i++)
		{
			MPI_Send(&vec[vecsize / ProcNum*i], vecsize / ProcNum, MPI_INT, i, 0, MPI_COMM_WORLD);
			cout << "send-" << ProcRank<<"->"<<i << endl;
		}
	}
	else
	{
		MPI_Recv(&vec[vecsize / (ProcNum*ProcRank)], vecsize / ProcNum, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		cout << "recv-" << ProcRank << "<-0" << endl;
		
	}
	_max = vec[0];
	cout  <<"working-"<< ProcRank << endl;
	/*if ((vecsize%ProcNum)&&(ProcRank==ProcNum-1))
	{
		for (int i = (vecsize / ProcNum)*ProcRank; i < vecsize; i++)
		{
			if (vec[i] > _max) _max = vec[i];
		}
	}
	for (int i = (vecsize/ProcNum)*ProcRank; i < (vecsize / ProcNum)*(ProcRank+1); i++)
	{
		if (vec[i] > _max) _max = vec[i];
	}*/
	_max = FindMax(vec, vecsize, ProcNum, ProcRank);
	MPI_Reduce(&_max, &max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	time = MPI_Wtime() - time;
	MPI_Finalize();
	cout <<"max="<< max << endl;
	//cout <<ProcRank << "_max=" << _max << endl;
	cout << "time=" << time << endl;
	return 0;
}