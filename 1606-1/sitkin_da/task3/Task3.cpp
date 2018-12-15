#include <mpi.h>
#include <cstdlib>
#include <time.h>
#include <Windows.h>
#include <queue>
#include <iostream>
#include <string>

#define INFINITI 10000000
#define WEIGHT 5
using namespace std;

void print_d(int* d, int size) 
{
	for (int i = 0; i < size; i++) 
	{
		cout << d[i] << " ";
	}
	cout << endl;
}
void print_graph(int** G, int size)
{
	for (int i = 0; i < size; i++) 
	{
		for (int j = 0; j < size; j++)
		{
			cout << G[i][j] << " ";
		}
		cout << endl;
	}
}

int* init_graph(int countEdge, int countVertex)
{
	int rank, procNum;

	MPI_Comm_size(MPI_COMM_WORLD, &procNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int* one_graph = new int[countVertex * countVertex];
	int** graph = new int*[countVertex];

	for (int i = 0; i < countVertex; i++) 
	{
		graph[i] = new int[countVertex];
	}
	srand(time(NULL));

	for (int i = 0; i < countVertex; i++) 
	{
		for (int j = 0; j < countVertex; j++)
		{
			if (i == j) 
			{
				graph[i][j] = 0;
			}
			else 
			{
				graph[i][j] = rand() % 3;
				if (graph[i][j] == 0) 
				{
					graph[i][j] = INFINITI;
				}
			}
			//cout << G[i][j] << " ";
		}
		//cout << endl;
	}
	for (int i = 0, t = 0; i < countVertex; i++) 
	{
		for (int j = 0; j < countVertex; j++)
		{
			one_graph[t] = graph[j][i];
			t++;
		}
	}
	return one_graph;
}

int* init_d(int size) 
{
	int* d = new int[size];
	for (int i = 0; i < size; i++) 
	{
		d[i] = INFINITI;
	}
	return d;
}

int* dijkstra(int* graph, int start, int count_vertex) 
{
	int* d = init_d(count_vertex);
	d[start] = 0;
	priority_queue<pair<int, int>>  queue;
	queue.push(make_pair(0, start));
	while (!queue.empty()) 
	{
		int v = queue.top().second, cur_d = queue.top().first;
		queue.pop();

		if (cur_d > d[v])  continue;

		for (int i = 0; i < count_vertex; ++i) 
		{
			int to = i, len = graph[i * count_vertex + v];
			if (d[v] + len < d[to]) 
			{
				d[to] = d[v] + len;
				queue.push(make_pair(d[to], to));
			}
		}
	}
	return d;
}

int* parallel_dijkstra(int* graph, int start, int count_vertex) 
{
	int* d = nullptr;
	int rank = 0, proc_num = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int* dist_d = new int[proc_num];
	int* dist_graph = new int[proc_num];
	int* count_element_d = new int[proc_num];
	int* count_element_graph = new int[proc_num];
	int part_size = count_vertex / (proc_num - 1);
	count_element_d[0] = 0; 
	dist_d[0] = 0;
	count_element_graph[0] = 0;
	dist_graph[0] = 0;
	for (int i = 1; i < proc_num - 1; i++)
	{
		count_element_d[i] = part_size;
		dist_d[i] = (i - 1) * part_size;
		count_element_graph[i] = count_vertex * part_size;
		dist_graph[i] = (i - 1) * count_vertex * part_size;
	}
	dist_d[proc_num - 1] = (proc_num - 2) * part_size;
	count_element_d[proc_num - 1] = count_vertex - (proc_num - 2) * part_size;
	dist_graph[proc_num - 1] = (proc_num - 2) * (count_vertex * part_size);
	count_element_graph[proc_num - 1] = count_vertex * count_vertex - (proc_num - 2) * (count_vertex*part_size);
	int* part_graph = new int[count_element_graph[rank]];
	int* part_d = new int[count_element_d[rank]];
	int flag = 1, flag_finish_find_current_min_distation_in_sigment_of_rank = 1;
	pair<int, int> current_vertex;
	pair<int, int> current_vertex_with_current_min_destation;
	MPI_Status st;
	if (rank == 0) 
	{
		MPI_Scatterv(graph, count_element_graph, dist_graph, MPI_INT, part_graph, count_element_graph[rank], MPI_INT, 0, MPI_COMM_WORLD); 
		d = init_d(count_vertex);
		d[start] = 0;
		priority_queue<pair<int, int>>  queue;
		queue.push(make_pair(0, start));

		int flag_rank_0 = 1;
		while (!queue.empty())
		{
			int v = queue.top().second;
			int	dest = queue.top().first;
			queue.pop();
			if (dest > d[v])
				continue;
			MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
			current_vertex = make_pair(dest, v);
			MPI_Bcast(&current_vertex, 1, MPI_2INT, 0, MPI_COMM_WORLD); 
			MPI_Scatterv(d, count_element_d, dist_d, MPI_INT, part_d, count_element_d[rank], MPI_INT, 0, MPI_COMM_WORLD);
			for (int i = 1; i < proc_num; i++)
			{
				flag_rank_0 = 1;
				while (flag_rank_0 != 0)
				{
					//std::cout << "wait" << std::endl;
					MPI_Recv(&flag_finish_find_current_min_distation_in_sigment_of_rank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &st);
					//std::cout << "get" << std::endl;
					flag_rank_0 = flag_finish_find_current_min_distation_in_sigment_of_rank;
					if (flag_rank_0 == 0) 
					{
						continue;
					}
					MPI_Recv(&current_vertex_with_current_min_destation, 1, MPI_2INT, i, 0, MPI_COMM_WORLD, &st);
					d[current_vertex_with_current_min_destation.second] = current_vertex_with_current_min_destation.first;

					queue.push(make_pair(current_vertex_with_current_min_destation.first, current_vertex_with_current_min_destation.second));
				}
			}
			flag_rank_0 = 1;
		}
		flag = 0;
		MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	else
	{
		MPI_Scatterv(graph, count_element_graph, dist_graph, MPI_INT, part_graph, count_element_graph[rank], MPI_INT, 0, MPI_COMM_WORLD);

		while (flag)
		{
			//std::cout << "succ" <<rank<< std::endl;
			MPI_Bcast(&flag, 1, MPI_INT, 0, MPI_COMM_WORLD);
			if (flag == 0) 
			{
				continue;
			}
			MPI_Bcast(&current_vertex, 1, MPI_2INT, 0, MPI_COMM_WORLD);
			MPI_Scatterv(d, count_element_d, dist_d, MPI_INT, part_d, count_element_d[rank], MPI_INT, 0, MPI_COMM_WORLD);

			for (int j = 0; j < count_element_d[rank]; j++) 
			{

				int to = j, len = part_graph[j * count_vertex + current_vertex.second];

				if (current_vertex.first + len < part_d[to])
				{
					flag_finish_find_current_min_distation_in_sigment_of_rank = 1;
					//std::cout << "send" << std::endl;
					MPI_Send(&flag_finish_find_current_min_distation_in_sigment_of_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
					current_vertex_with_current_min_destation.first = current_vertex.first + len;
					current_vertex_with_current_min_destation.second = to + dist_d[rank];

					MPI_Send(&current_vertex_with_current_min_destation, 1, MPI_2INT, 0, 0, MPI_COMM_WORLD);
				}
			}

			flag_finish_find_current_min_distation_in_sigment_of_rank = 0;
			MPI_Send(&flag_finish_find_current_min_distation_in_sigment_of_rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}
	return d;
}
int main(int argc, char *argv[])
{
	int proc_num, proc_rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

	const int count_vertex = stoi(std::string(argv[1]));
	srand(time(NULL));
	const int count_edge = (count_vertex - 1) + rand() % ((count_vertex * (count_vertex - 1)) / 2);

	int* d_step = nullptr;
	int* d_parallel = nullptr;
	int* graph = nullptr;
	double t1_step = 0.0, t2_step = 0.0, t1_parallel = 0.0, t2_parallel = 0.0;
	int start = rand() % (count_vertex - 1);
	if (proc_rank == 0) 
	{
		graph = init_graph(count_edge, count_vertex);
	}
	if (proc_rank == 0)
	{
		t1_parallel = MPI_Wtime();
	}
	d_parallel = parallel_dijkstra(graph, start, count_vertex);
	MPI_Barrier(MPI_COMM_WORLD);
	if (proc_rank == 0) 
	{
		t2_parallel = MPI_Wtime();
		cout << endl << " Parallel Dijkstra time: " << t2_parallel - t1_parallel << endl << endl;
		t1_step = MPI_Wtime();
		d_step = dijkstra(graph, start, count_vertex);
		t2_step = MPI_Wtime();
		cout << endl << " Dijkstra time: " << t2_step - t1_step << endl << endl;
	}
	MPI_Finalize();
}