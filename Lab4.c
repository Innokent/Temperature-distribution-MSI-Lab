/*
3. Разработать средствами MPI параллельную программу решения двухмерной нестационарной краевой задачи методом конечных разностей с использованием явной вычислительной схемы. Объект моделирования - прямоугольная пластина постоянной толщины. Подробности постановки подобной задачи даны ниже. Возможны граничные условия первого и второго рода в различных узлах расчетной сетки. Временной интервал моделирования и количество (кратное 8) узлов расчетной сетки - параметры программы. Программа должна демонстрировать ускорение по сравнению с последовательным вариантом. Предусмотреть визуализацию результатов посредством утилиты gnuplot.
*/ 
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>

#define k 1 // Коэффицент теплопроводности
#define Gt 0 // Приведенная скорость взаимного превращения тепловой энергии в другие виды энергии
#define dx 2 // Шаг по Х
#define dy 2 // Шаг по Y
#define L 100 // Длинна стержня
#define S 10 // Площадь поперечного сечения стержня
#define boundaryConditionFirstTop 10
#define boundaryConditionSecondTop 0.1
#define boundaryConditionFirstBot 20
#define boundaryConditionSecondBot 0.2
#define boundaryConditionFirstLeft 30
#define boundaryConditionSecondLeft 0.3
#define boundaryConditionFirstRight 40
#define boundaryConditionSecondRight 0.4
#define maxTime 500
#define currMode 1


typedef struct 
{
	pthread_t tid;
	int firstIndexStart;
	int firstIndexEnd;
	int secondIndexStart;
	int secondIndexEnd;
	double dt;
	int n;
	int m;
} Thread_param;

pthread_mutex_t mutx1;
pthread_barrier_t barr1;
Thread_param *threads;
double *prevLayer;
double *currLayer;

double calcNext(double T01, double T11, double T21, double T10, double T12, double dt)
{
	double numerator1 = T01 - 2 * T11 + T21;
	double numerator2 = T10 - 2 * T11 + T12;
	double fraction1 = numerator1 / (dx * dx);
	double fraction2 = numerator2 / (dy * dy);
	return (dt * (k * (fraction1 + fraction2) + Gt) + T11);
}

double calcBoundary(int i, int j, int n, int m, double dt, double T, int mode)
{
	if(mode == 1) // Граничные условия первого рода
	{
		
		if(i == 0)
		{
			return boundaryConditionFirstTop;
		}
		if(i == n - 1)
		{
			return boundaryConditionFirstBot;
		}
		if(j == 0)
		{
			return boundaryConditionFirstLeft;
		}
		if(j == m - 1)
		{
			return boundaryConditionFirstRight;
		}
	}
	else // Граничные условия второго рода
	{
		if(i == 0)
		{
			return (boundaryConditionSecondTop * dt / dx + 1) * T;
		}
		if(i == n - 1)
		{
			return (boundaryConditionSecondBot * dt / dx + 1) * T;
		}
		if(j == 0)
		{
			return (boundaryConditionSecondLeft * dt / dy + 1) * T;
		}
		if(j == m - 1)
		{
			return (boundaryConditionSecondRight * dt / dy + 1) * T;
		}
	}
}

void *solver(void *arg_p)
{
	Thread_param *params = (Thread_param *) arg_p;
	double T01, T11, T21, T10, T12;
	int K = 0;
	double *temp = prevLayer;
	
	for(double t = 0.0 + params->dt; t <= maxTime; t += params->dt)
	{
		pthread_barrier_wait(&barr1);
		for(int i = params->firstIndexStart; i <= params->firstIndexEnd; i++ )
		{
			for(int j = params->secondIndexStart; j <= params->secondIndexEnd; j++)
			{
				if((i != 0) && (i != (params->n - 1)) && (j != 0) && (j != (params->m - 1)))
				{
					T01 = prevLayer[params->n * j + i - 1];
					T11 = prevLayer[params->n * j + i];
					T21 = prevLayer[params->n * j + i + 1];
					T10 = prevLayer[params->n * (j - 1) + i];
					T12 = prevLayer[params->n * (j + 1) + i];
					currLayer[params->n * j + i] = calcNext(T01, T11, T21, T10, T12, params->dt);
				}
				else
				{
					currLayer[params->n * j + i] = 
							calcBoundary(i, j, params->n, params->m, params->dt, prevLayer[params->n * j + i], currMode);
				}
			}
		}
		K++;
		pthread_barrier_wait(&barr1);
		
		if(!pthread_mutex_trylock(&mutx1))
		{
			temp = prevLayer;
			prevLayer = currLayer;
			currLayer = temp;
			
			FILE *output = fopen("output.txt", "a");
			for(int i = 0; i < params->n; i++)
			{
				for(int j = 0; j < params->m; j++)
				{
					fprintf(output, "%d %d %lf\n", i, j, currLayer[params->n * j + i]);
				}
			}
			fprintf(output, "\n\n");
			fclose(output);

			pthread_mutex_unlock(&mutx1);
		}
	}	
}

int main(int argc, char* argv[])
{
	pthread_attr_t pattr;
	
	if(argc != 5) // Проверяем количество параметров
	{
		puts("Ошибка: Неверное количесвто параметров");
		exit(1);
	}
	
	int c_threads = atoi(argv[1]);
	double dt = atof(argv[2]);
	int N = atoi(argv[3]);
	int M = atoi(argv[4]);
	int K = (int)(maxTime / dt) + 1;
	
	if((N * M) % 8) // Проверям кратность числа узлов восьми
	{
		puts("Ошибка: Количество узлов сетки должно быть кратно 8");
		exit(3);
	}
	
	if(N % c_threads) // Проверяем условие делимости нацело числа узлов и числа потоков
	{
		puts("Ошибка: Количество узов сетки по одной координте не делится нацело на количество потоков");
		exit(4);
	}

	N += 2; // Учитываем граничные условия для верха и низа пластины
	M += 2; // Учитываем граничные условия для левого и правого края пластины

	prevLayer = (double *) calloc(N * M, sizeof(double)); // Выделяем память под хранение предыдущего слоя
	currLayer = (double *) calloc(N * M, sizeof(double)); // Выделяем память под хранение текущего слоя

	if((!prevLayer) || (!currLayer)) // Проверяем удалось ли выделить память под массив
	{
		puts("Ошибка: Не удалось выделить память");
		exit(2);
	}
	
	for(int i = 0; i < N; i++) // Заполняем узлы балки значениями в начальный момент времни (t = 0)
	{
		for(int j = 0; j < M; j++)
		{
			if((i != 0) && (i != N - 1) && (j != 0) && (j != M - 1))
			{
				prevLayer[N * j + i] = 0.0;
			}
			else
			{
				prevLayer[N * j + i] = calcBoundary(i, j, N, M, dt, 0, 1);
			}
		}
	}
	
	FILE *output = fopen("output.txt", "w");
	for(int i = 0; i < N; i++)
	{
		for(int j = 0; j < M; j++)
		{
			fprintf(output, "%d %d %lf\n", i, j, prevLayer[N * j + i]);
		}
	}
	fprintf(output, "\n\n");
	fclose(output);

	pthread_attr_init(&pattr);
	pthread_attr_setscope(&pattr, PTHREAD_SCOPE_SYSTEM);
	pthread_attr_setdetachstate(&pattr, PTHREAD_CREATE_JOINABLE);

	pthread_mutex_init(&mutx1, NULL);
	pthread_barrier_init(&barr1, NULL, c_threads);

	threads = (Thread_param *) calloc(c_threads, sizeof(Thread_param)); // Выделяем память под потоки

	for(int i = 0; i < c_threads; i++) // Инициализируем атрибуты потоков
	{
		threads[i].n = N;
		threads[i].m = M;
		threads[i].dt = dt;
		threads[i].firstIndexStart = 1 + i * ((N - 2) / c_threads);
		if(i == 0)
		{
			threads[i].firstIndexStart--;
		}
		threads[i].firstIndexEnd = (i + 1) * ((N - 2) / c_threads);
		if(i == c_threads - 1)
		{
			threads[i].firstIndexEnd++;
		}
		threads[i].secondIndexStart = 0;
		threads[i].secondIndexEnd = M - 1;
		
	}
	
	for(int i = 0; i < c_threads; i++)
	{
		if(pthread_create(&threads[i].tid, &pattr, solver, (void *) &(threads[i]))) //Создаем потоки
		{
			puts("Ошибка: Не удалось создать поток");
			exit(5);
		}
	}

	pthread_join(threads[c_threads - 1].tid, NULL);

	FILE *gnuplot = popen("gnuplot -persist", "w");

	if(gnuplot == NULL)
	{
		puts("Ошибка: Не удалось запустить утилиту Gnuplot");
		exit(6);
	}

	fputs("set dgrid3d\n", gnuplot);
	fputs("set hidden3d\n", gnuplot);

	for(int i = 0; i < K - 1; i++)
	{
		fprintf(gnuplot, "splot 'output.txt' index %d with lines\n", i);
		usleep(10000);
		fflush(gnuplot);
	}
	
	pclose(gnuplot);
	pthread_mutex_destroy(&mutx1);
	pthread_barrier_destroy(&barr1);
	free(prevLayer);
	free(currLayer);
	free(threads);

	return 1;
}
