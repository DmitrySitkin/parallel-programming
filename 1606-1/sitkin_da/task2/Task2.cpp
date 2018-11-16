#include"mpi.h"
#include<iostream>
void my_bcast(void* data, int count, MPI_Datatype datatype, int root,
	MPI_Comm communicator) {
	int world_rank;
	MPI_Comm_rank(communicator, &world_rank);
	int world_size;
	MPI_Comm_size(communicator, &world_size);

	if (world_rank == root) {
		int i;
		for (i = 0; i < world_size; i++) {
			if (i != world_rank) {
				MPI_Send(data, count, datatype, i, 0, communicator);
			}
		}
	}
	else {
		MPI_Recv(data, count, datatype, root, 0, communicator,
			MPI_STATUS_IGNORE);
	}
}
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
	int _max = 0;
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
	int n, ProcNum, ProcRank, *vec, max, _max;
	int vecsize = 99900000;
	double time = 0.0;
	MPI_Status status;
	srand(1);
	vec = CreateRandVec(vecsize);
	vec[9999] = 999999999;
	MPI_Init(&argc, &argv);
	time = MPI_Wtime();
	//cout << "2" << endl;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	std::cout << "proc" << ProcRank << std::endl;
	my_bcast(&vecsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	std::cout << "Bcast-" << ProcRank << std::endl;
	if (ProcRank == 0)
	{
		for (int i = 1; i < ProcNum; i++)
		{
			MPI_Send(&vec[vecsize / ProcNum*i], vecsize / ProcNum, MPI_INT, i, 0, MPI_COMM_WORLD);
			std::cout << "send-" << ProcRank << "->" << i << std::endl;
		}
	}
	else
	{
		MPI_Recv(&vec[vecsize / (ProcNum*ProcRank)], vecsize / ProcNum, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		std::cout << "recv-" << ProcRank << "<-0" << std::endl;

	}
	_max = vec[0];
	std::cout << "working-" << ProcRank << std::endl;
	_max = FindMax(vec, vecsize, ProcNum, ProcRank);
	MPI_Reduce(&_max, &max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	time = MPI_Wtime() - time;
	MPI_Finalize();
	std::cout << "max=" << max << std::endl;
	//cout <<ProcRank << "_max=" << _max << endl;
	std::cout << "time=" << time << std::endl;
	return 0;
}
