#include <stdio.h>
#include <omp.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include <x86intrin.h>

#include "sorting.h"

/* 
   Merge two sorted chunks of array T!
   The two chunks are of size size
   First chunck starts at T[0], second chunck starts at T[size]
*/
void merge(uint64_t *T, const uint64_t size)
{
	uint64_t *X = (uint64_t *)malloc(2 * size * sizeof(uint64_t));

	uint64_t i = 0;
	uint64_t j = size;
	uint64_t k = 0;

	while ((i < size) && (j < 2 * size))
	{
		if (T[i] < T[j])
		{
			X[k] = T[i];
			i = i + 1;
		}
		else
		{
			X[k] = T[j];
			j = j + 1;
		}
		k = k + 1;
	}

	if (i < size)
	{
		for (; i < size; i++, k++)
		{
			X[k] = T[i];
		}
	}
	else
	{
		for (; j < 2 * size; j++, k++)
		{
			X[k] = T[j];
		}
	}

	memcpy(T, X, 2 * size * sizeof(uint64_t));
	free(X);

	return;
}

/* 
   merge sort -- sequential, parallel -- 
*/

void sequential_merge_sort(uint64_t *T, const uint64_t size)
{
	if (size < 2)
		return;
	uint64_t subsize = size / 2;
	sequential_merge_sort(T, subsize);
	sequential_merge_sort(T + subsize, subsize);
	merge(T, subsize);

	return;
}

void parallel_merge_sort(uint64_t *T, const uint64_t size)
{
	if (size <= 64){
		sequential_merge_sort(T,size);
		return;
	}

	#pragma omp parallel
	{
		#pragma omp single
		{
			uint64_t subsize = size / 2;
		#pragma omp task depend(inout: T)
			{parallel_merge_sort(T, subsize);}
		#pragma omp task depend(inout: T)
			{parallel_merge_sort(T + subsize, subsize);}
		#pragma omp taskwait
			merge(T, subsize);
		}
	}
	
}


void parallel_merge_sort_iter(uint64_t *T, const uint64_t size)
{
	int threashold = 64;
	if (size <= threashold){
		sequential_merge_sort(T,size);
		return;
	}
	#pragma omp parallel
	{
		#pragma omp single
		{
			for(int i=threashold; i<=size; i*=2){
				#pragma omp taskwait
				for(int j=0; j<size; j+=i){
					#pragma omp task depend(inout: T)
					{
						if(i<=threashold){
							sequential_merge_sort(T+j,i);
						}else{
							merge(T+j,i/2);
						}
					}
				}
			}
		}
	}
}


int main(int argc, char **argv)
{
	struct timespec begin, end;
	double seconds;
	double nanoseconds;

	unsigned int exp;

	/* the program takes one parameter N which is the size of the array to
       be sorted. The array will have size 2^N */
	if (argc != 2)
	{
		fprintf(stderr, "merge.run N \n");
		exit(-1);
	}

	uint64_t N = 1 << (atoi(argv[1]));
	//printf("%d",atoi(argv[1]));
	/* the array to be sorted */
	uint64_t *X = (uint64_t *)malloc(N * sizeof(uint64_t));

	printf("--> Sorting an array of size %lu\n", N);
#ifdef RINIT
	printf("--> The array is initialized randomly\n");
#endif

	for (exp = 0; exp < NBEXPERIMENTS; exp++)
	{
#ifdef RINIT
		init_array_random(X, N);
#else
		init_array_sequence(X, N);
#endif

		clock_gettime(CLOCK_MONOTONIC, &begin);

		sequential_merge_sort(X, N);

		clock_gettime(CLOCK_MONOTONIC, &end);

		seconds = end.tv_sec - begin.tv_sec;
		nanoseconds = end.tv_nsec - begin.tv_nsec;

		experiments[exp] = seconds + nanoseconds * 1e-9;

		// verifying that X is properly sorted //
#ifdef RINIT
		if (!is_sorted(X, N))
		{
			print_array(X, N);
			fprintf(stderr, "ERROR: the sequential sorting of the array failed\n");
			exit(-1);
		}
#else
		if (!is_sorted_sequence(X, N))
		{
			print_array(X, N);
			fprintf(stderr, "ERROR: the sequential sorting of the array failed\n");
			exit(-1);
		}
#endif
	}

	//printf("\n mergesort serial \t\t\t %.6lf seconds\n\n", average_time());
	printf(";%.3lf",1000*average_time());

	for (exp = 0; exp < NBEXPERIMENTS; exp++)
	{
#ifdef RINIT
		init_array_random(X, N);
#else
		init_array_sequence(X, N);
#endif

		clock_gettime(CLOCK_MONOTONIC, &begin);

		parallel_merge_sort(X, N);

		clock_gettime(CLOCK_MONOTONIC, &end);

		seconds = end.tv_sec - begin.tv_sec;
		nanoseconds = end.tv_nsec - begin.tv_nsec;

		experiments[exp] = seconds + nanoseconds * 1e-9;

		/* verifying that X is properly sorted */
#ifdef RINIT
		if (!is_sorted(X, N))
		{
			print_array(X, N);
			fprintf(stderr, "ERROR: the parallel sorting of the array failed\n");
			exit(-1);
		}
#else
		if (!is_sorted_sequence(X, N))
		{
			print_array(X, N);
			fprintf(stderr, "ERROR: the parallel sorting of the array failed\n");
			exit(-1);
		}
#endif
	}

	printf("\n mergesort parallel \t\t\t %.6lf seconds\n\n", average_time());
	//printf(";%.3lf",1000*average_time());


	for (exp = 0; exp < NBEXPERIMENTS; exp++)
	{
#ifdef RINIT
		init_array_random(X, N);
#else
		init_array_sequence(X, N);
#endif

		clock_gettime(CLOCK_MONOTONIC, &begin);

		parallel_merge_sort_iter(X, N);

		clock_gettime(CLOCK_MONOTONIC, &end);

		seconds = end.tv_sec - begin.tv_sec;
		nanoseconds = end.tv_nsec - begin.tv_nsec;

		experiments[exp] = seconds + nanoseconds * 1e-9;

		/* verifying that X is properly sorted */
#ifdef RINIT
		if (!is_sorted(X, N))
		{
			print_array(X, N);
			fprintf(stderr, "ERROR: the parallel sorting of the array failed\n");
			exit(-1);
		}
#else
		if (!is_sorted_sequence(X, N))
		{
			print_array(X, N);
			fprintf(stderr, "ERROR: the parallel sorting of the array failed\n");
			exit(-1);
		}
#endif
	}

	printf("\n mergesort parallel iterative \t\t\t %.6lf seconds\n\n", average_time());
	//printf(";%.3lf",1000*average_time());


	/* print_array (X, N) ; */

	/* before terminating, we run one extra test of the algorithm */
	uint64_t *Y = (uint64_t *)malloc(N * sizeof(uint64_t));
	uint64_t *Z = (uint64_t *)malloc(N * sizeof(uint64_t));

#ifdef RINIT
	init_array_random(Y, N);
#else
	init_array_sequence(Y, N);
#endif

	memcpy(Z, Y, N * sizeof(uint64_t));

	sequential_merge_sort(Y, N);
	parallel_merge_sort(Z, N);

	if (!are_vector_equals(Y, Z, N))
	{
		fprintf(stderr, "ERROR: sorting with the sequential and the parallel algorithm does not give the same result\n");
		exit(-1);
	}

	free(X);
	free(Y);
	free(Z);
}
