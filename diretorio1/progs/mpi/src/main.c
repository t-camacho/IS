/*
 * Copyright(C) 2014 Pedro H. Penna <pedrohenriquepenna@gmail.com>
 *
 * Integer-Sort Benchmark Kernel.
 */

#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <timer.h>
#include <util.h>
#include "is.h"

#define NUM_BUCKETS 8192
#define SZ_DARRAY_TYPE 2
/*
 * Problem.
 */
struct problem
{
	int n;   /* Number of elements. */
	int max; /* Maximum number.     */
};

/* Problem sizes. */
static struct problem tiny     = {   8388608, (1 << 20) };
static struct problem small    = {  16777216, (1 << 20) };
static struct problem standard = {  33554432, (1 << 20) };
static struct problem large    = {  67108864, (1 << 20) };
static struct problem huge     = { 134217728, (1 << 20) };

/* Be verbose? */
int verbose = 0;

/* Number of threads. */
int nthreads = 1;

/* Seed number. */
static int seed = 0;

/* Problem. */
static struct problem *p = &tiny;

MPI_Datatype MPI_DARRAY;
/*
 * Prints program usage and exits.
 */
static void usage(void)
{
	printf("Usage: mpirun -n <value> ./bin/is.intel [options]\n");
	printf("Options:\n");
	printf("  --help             Display this information and exit\n");
	printf("  --class <name>     Set problem class:\n");
	printf("                       - small\n");
	printf("                       - standard\n");
	printf("                       - large\n");
	printf("                       - huge\n");
	printf("  --verbose          Be verbose\n");
	exit(0);
}

void write_file(int *array, char *n) {
	FILE *file = fopen(n, "w+");

	for (int i = 0; i < 100; i++)
	{
		fprintf(file, "%i ", array[i]);
	}

	fclose(file);
}

void build_mpi_type()
{
	int blocklens[SZ_DARRAY_TYPE];
	MPI_Aint indices[SZ_DARRAY_TYPE];
	MPI_Datatype old_types[SZ_DARRAY_TYPE];

	for(int i = 0; i < SZ_DARRAY_TYPE; i++)
		old_types[i] = MPI_INT;

	for(int i = 0; i < SZ_DARRAY_TYPE; i++)
		blocklens[i] = 1;

	indices[0] = 0;
	indices[1] = 4;

	MPI_Type_create_struct(SZ_DARRAY_TYPE, blocklens, indices, old_types, &MPI_DARRAY);
	MPI_Type_commit(&MPI_DARRAY);
}

int get_max(int start, int end, int* array) {
	int max = INT_MIN;

	for (int i = start; i <= end; i++)
	{
		/* Found. */
		if (array[i] > max)
			max = array[i];
	}

	return max;
}

/*
 * Reads command line arguments.
 */
static void readargs(int argc, char **argv)
{
	int i;          /* Loop index.         */
	char *arg; /* Working argument. */
	int state; /* Processing state. */

	/* State values. */
	#define READ_ARG     0 /* Read argument.         */
	#define SET_NTHREADS 1 /* Set number of threads. */
	#define SET_CLASS    2 /* Set problem class.     */

	state = READ_ARG;

	/* Read command line arguments. */
	for (i = 1; i < argc; i++)
	{
		arg = argv[i];

		/* Set value. */
		if (state != READ_ARG)
		{
			switch (state)
			{
				/* Set problem class. */
				case SET_CLASS :
					if (!strcmp(argv[i], "tiny"))
						p = &tiny;
					else if (!strcmp(argv[i], "small"))
						p = &small;
					else if (!strcmp(argv[i], "standard"))
						p = &standard;
					else if (!strcmp(argv[i], "large"))
						p = &large;
					else if (!strcmp(argv[i], "huge"))
						p = &huge;
					else
						usage();
					state = READ_ARG;
					break;

				/* Set number of threads. */
				case SET_NTHREADS :
					nthreads = atoi(arg);
					state = READ_ARG;
					break;

				default:
					usage();
			}
			continue;
		}

		/* Parse argument. */
		if (!strcmp(arg, "--verbose"))
			verbose = 1;
		else if (!strcmp(arg, "--nthreads"))
			state = SET_NTHREADS;
		else if (!strcmp(arg, "--class"))
			state = SET_CLASS;
		else
			usage();
	}

	/* Invalid argument(s). */
	if (nthreads < 1)
		usage();
}

/*
 * Runs benchmark.
 */
int main(int argc, char **argv)
{
	int i, j, k;          /* Loop index.         */
	int *array;           /* Array to be sorted. */
	int my_rank, comm_sz; /* mpi                 */
	int sz;               /* size subarray       */
	int start, end;       /* start/end subarray  */
	struct darray **buckets; /* Buckets.            */

#ifdef _XEON_PHI_
	double power;
#endif

	readargs(argc, argv);

	MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	sz = (p->n)/comm_sz;
	start = my_rank * sz;
	end   = start + sz - 1;

	array = smalloc(p->n*sizeof(int));

	build_mpi_type();

	if(my_rank == 0)
	{
		double num;       /* Normal number.      */
		uint64_t t_end;   /* End time.           */
		uint64_t t_start; /* Start time.         */
		int max, range;

		int *indexes;            /* Index for buckets.   */

		timer_init();
		srandnum(seed);

		// write_file(array, "output/unordered.out");

		/* Benchmark initialization. */
		if (verbose)
			printf("initializing...\n");

		t_start = timer_get();
		for (i = 0; i < p->n; i++)
		{
			num = normalnum(0, (p->max >> 4));
			array[i] = (int)((num < 0) ? -num : num) + 1;
		}
		t_end = timer_get();

		if (verbose)
			printf("  time spent: %f\n", timer_diff(t_start, t_end)*MICROSEC);

		#ifdef _XEON_PHI_
			power_init();
		#endif

		/* Cluster data. */
		if (verbose)
			printf("sorting...\n");

		t_start = timer_get();

		indexes = smalloc(NUM_BUCKETS*sizeof(int));

		/* Create buckets. */
		buckets = smalloc(NUM_BUCKETS*sizeof(struct darray *));
		for (i = 0; i < NUM_BUCKETS; i++)
			buckets[i] = darray_create(p->n/NUM_BUCKETS);

		/*max*/
		for(int i = 1; i < comm_sz; i++) {
      MPI_Send(array, p->n, MPI_INT, i, 0, MPI_COMM_WORLD);
    }

		max = get_max(start, end, array);

		for(int i = 1; i < comm_sz; i++) {
			int buff;
      MPI_Recv(&buff, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (buff > max)
				max = buff;
		}

		/*bucket*/
		range = max/NUM_BUCKETS;

		for (i = 0; i < p->n; i++)
		{
			j = array[i]/range;
			if (j >= NUM_BUCKETS)
				j = NUM_BUCKETS - 1;
			darray_append(buckets[j], array[i]);
		}

		/* Sort Each bucket. */
		sz = NUM_BUCKETS/(comm_sz-1);
		for(i = 1; i < comm_sz; i++) {
			start = (i-1) * sz;
			end   = start + sz - 1;
			for (int k = start; k <= end; k++)
			{
				MPI_Send(buckets[k], 1, MPI_DARRAY, i, 0, MPI_COMM_WORLD);
				MPI_Send(buckets[k]->elements, buckets[k]->size, MPI_INT, i, 0, MPI_COMM_WORLD);
			}
		}

		/* Destroy buckets. */
		for (i = 0; i < NUM_BUCKETS; i++)
			darray_destroy(buckets[i]);

		/* Create buckets. */
		buckets = smalloc(NUM_BUCKETS*sizeof(struct darray *));
		for (i = 0; i < NUM_BUCKETS; i++)
			buckets[i] = darray_create(p->n/NUM_BUCKETS);

		k = 0;
		for(i = 1; i < comm_sz; i++) {
			for(j = 0; j < sz; j++) {
				MPI_Recv(buckets[k], 1, MPI_DARRAY, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				buckets[k]->elements = (int*) malloc(buckets[k]->size*sizeof(int));
				MPI_Recv(buckets[k]->elements, buckets[k]->size, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				k++;
			}
		}

		/* Build indexes. */
		indexes[0] = 0;
		for (i = 1; i < NUM_BUCKETS; i++)
			indexes[i] = indexes[i - 1] + darray_size(buckets[i]);

		/* Rebuild array. */
		for (i = 0; i < NUM_BUCKETS; i++)
		{
			int k = indexes[i];

			for (j = 0; j < darray_size(buckets[i]); j++)
				array[k + j] = darray_get(buckets[i], j);
		}

		/* House keeping. */
		for (i = 0; i < NUM_BUCKETS; i++)
			darray_destroy(buckets[i]);

		t_end = timer_get();

		printf("timing statistics:\n");
		printf("  total time:    %f\n",  timer_diff(t_start, t_end)*MICROSEC);

		// write_file(array, "output/ordered.out");

		#ifdef _XEON_PHI_
			power = power_end();
		#endif

		#ifdef _XEON_PHI_
			printf("  average power: %f\n", power*0.000001);
		#endif

	}else {
		MPI_Recv(array, p->n, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		int _max = get_max(start, end, array);

		MPI_Send(&_max, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

		/* Create buckets. */
		buckets = smalloc(NUM_BUCKETS/(comm_sz-1)*sizeof(struct darray *));
		for (i = 0; i < NUM_BUCKETS; i++)
			buckets[i] = darray_create(p->n/NUM_BUCKETS);

		for(i = 0; i < NUM_BUCKETS/(comm_sz-1); i++)
		{
			MPI_Recv(buckets[i], 1, MPI_DARRAY, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			buckets[i]->elements = (int*) malloc(buckets[i]->size*sizeof(int));
			MPI_Recv(buckets[i]->elements, buckets[i]->size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (darray_size(buckets[i]) > 0)
					sort(buckets[i]);
		}

		for(i = 0; i < NUM_BUCKETS/(comm_sz-1); i++)
		{
			MPI_Send(buckets[i], 1, MPI_DARRAY, 0, 0, MPI_COMM_WORLD);
			MPI_Send(buckets[i]->elements, buckets[i]->size, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}

		/* House keeping. */
		for (i = 0; i < NUM_BUCKETS/(comm_sz-1); i++)
			darray_destroy(buckets[i]);
	}

	/* House keeping. */
	free(array);

	MPI_Finalize();

	return (0);
}
