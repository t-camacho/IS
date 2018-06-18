/*
 * Copyright(C) 2014 Pedro H. Penna <pedrohenriquepenna@gmail.com>
 *
 * Integer-Sort Benchmark Kernel.
 */

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <timer.h>
#include <util.h>
#include "is.h"

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

/* Seed number. */
static int seed = 0;

/* Problem. */
static struct problem *p = &tiny;

/*
 * Prints program usage and exits.
 */
static void usage(void)
{
	printf("Usage: ./bin/is.intel [options]\n");
	printf("Brief: Integer Sort Benchmark Kernel\n");
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

/*
 * Reads command line arguments.
 */
static void readargs(int argc, char **argv)
{
	int i;     /* Loop index.       */
	char *arg; /* Working argument. */
	int state; /* Processing state. */

	/* State values. */
	#define READ_ARG     0 /* Read argument.         */
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

				default:
					usage();
			}
			continue;
		}

		/* Parse argument. */
		if (!strcmp(arg, "--verbose"))
			verbose = 1;
		else if (!strcmp(arg, "--class"))
			state = SET_CLASS;
		else
			usage();
	}
}

void write_file(int *array, char *n) {
	FILE *file = fopen(n, "w+");

	for (int i = 0; i < 100; i++)
	{
		fprintf(file, "%i ", array[i]);
	}

	fclose(file);
}

/*
 * Runs benchmark.
 */
int main(int argc, char **argv)
{
	int i;          /* Loop index.         */
	int *a;         /* Array to be sorted. */
	double num;     /* Normal number.      */
	uint64_t end;   /* End time.           */
	uint64_t start; /* Start time.         */

#ifdef _XEON_PHI_
	double power;
#endif

	readargs(argc, argv);

	timer_init();
	srandnum(seed);

	/* Benchmark initialization. */
	if (verbose)
		printf("initializing...\n");

	start = timer_get();
	a = smalloc(p->n*sizeof(int));
	for (i = 0; i < p->n; i++)
	{
		num = normalnum(0, (p->max >> 4));
		a[i] = (int)((num < 0) ? -num : num) + 1;
	}
	end = timer_get();

	if (verbose)
		printf("  time spent: %f\n", timer_diff(start, end)*MICROSEC);

	// write_file(a, "output/unordered.out");

#ifdef _XEON_PHI_
	power_init();
#endif

	/* Cluster data. */
	if (verbose)
		printf("sorting...\n");

	start = timer_get();
	integer_sort(a, p->n);
	end = timer_get();

	printf("timing statistics:\n");
	printf("  total time:    %f\n",  timer_diff(start, end)*MICROSEC);

	// write_file(a, "output/ordered.out");

#ifdef _XEON_PHI_
	power = power_end();
#endif

#ifdef _XEON_PHI_
	printf("  average power: %f\n", power*0.000001);
#endif

	/* House keeping. */
	free(a);

	return (0);
}
