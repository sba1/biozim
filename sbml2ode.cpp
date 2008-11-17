#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Own headers */
#include "simulation.h"

/********************************************************/

static const char *model_filename;
static int force_interpreted;
static int stiff;
static int stochastic;
static int output_time;
static double maxtime;
static double error;
static int sample_steps;

int verbose; /* used by other modules */


/*************************************************
 Displays the program's usage and exits
*************************************************/
static void usage(char *name)
{
	fprintf(stderr, "usage: %s [OPTIONS] SBML-File ...\n"
			"Loads the given SBML-File and performs a simulation run.\n"
			"Specify '-' to read from stdin.\n"
			"\t-h, --help               show this help and quit.\n"
			"\t    --error              specifies the error (absolute and relative) given in double.\n"
			"\t    --force-interpreted  forces the interpreted calculation of the rhs.\n"
			"\t    --maxtime            specifies the end time (defaults to 1).\n"
			"\t    --output-time        outputs the current time (stderr).\n"
			"\t    --sample-steps       the number of sample steps (defaults to 5000).\n"
			"\t    --stiff              use solver for stiff ODEs.\n"
			"\t    --stochastic         apply stochastic simulation.\n"
			"\t    --verbose            verbose output.\n",
			name);

	exit(1);
}


/*************************************************
 Parse command line args
*************************************************/
static void parse_args(int argc, char *argv[])
{
	int i;
	int filename_given = 0;

	maxtime = 1.0;
	error = 1e-10;
	sample_steps = 5000;

	for (i=1;i<argc;i++)
	{
		if (!strcmp(argv[i],"-h") || !strcmp(argv[i],"--help"))
		{
			usage(argv[0]);
			exit(-1);
		} else if (!strcmp(argv[i],"--verbose"))
		{
			verbose = 1;
		} else if (!strcmp(argv[i],"--force-interpreted"))
		{
			force_interpreted = 1;
		} else if (!strcmp(argv[i],"--output-time"))
		{
			output_time = 1;
		} else if (!strcmp(argv[i],"--stiff"))
		{
			stiff = 1;
		} else if (!strcmp(argv[i],"--stochastic"))
		{
			stochastic = 1;
		} else if (!strcmp(argv[i], "--sample-steps"))
		{
			char *nr_arg;

			if (argv[i][9]=='=')
				nr_arg = &argv[i][10];
			else
			{
				nr_arg = argv[i+1];
				i++;

				if (i>=argc)
				{
					fprintf(stderr,"The --sample-steps option needs an integer argument.\n");
					exit(-1);
				}
			}
			sample_steps = strtod(nr_arg, NULL);
		} else if (!strcmp(argv[i],"--maxtime"))
		{
			char *nr_arg;

			if (argv[i][9]=='=')
				nr_arg = &argv[i][10];
			else
			{
				nr_arg = argv[i+1];
				i++;

				if (i>=argc)
				{
					fprintf(stderr,"The --maxtime option needs an integer argument.\n");
					exit(-1);
				}
			}
			maxtime = strtod(nr_arg, NULL);
		} else if (!strcmp(argv[i],"--error"))
		{
			char *nr_arg;

			if (argv[i][9]=='=')
				nr_arg = &argv[i][10];
			else
			{
				nr_arg = argv[i+1];
				i++;

				if (i>=argc)
				{
					fprintf(stderr,"The --error option needs an argument.\n");
					exit(-1);
				}
			}
			error = strtod(nr_arg, NULL);
		}

		else
		{
			filename_given = 1;
			if (strcmp(argv[i],"-"))
				model_filename = argv[i];
		}
	}

	if (!filename_given)
	{
		fprintf(stderr,"No filename has been specified!\n");
		usage(argv[0]);
		exit(-1);
	}

	if (sample_steps < 1) sample_steps = 1;
}

/**********************************************************
 Our sampling function.
***********************************************************/
int sample(double time, int num_values, double *values)
{
	if (output_time)
		fprintf(stderr,"%g\n",time);

	printf("%g",time);

	for (int i=0;i<num_values;i++)
	{
		printf("\t%.12g",values[i]);
	}

	return 1;
}

/**********************************************************
 Our sampling function for string. This is called always
 in sync with the above function.
***********************************************************/
int sample_strings(double time, int num_strings, char **values)
{
	for (int i=0;i<num_strings;i++)
	{
		printf("\t\"%s\"",values[i]);
	}
	printf("\n");

	return 1;
}

/**********************************************************
 Main Entry
***********************************************************/
int main(int argc, char **argv)
{
	struct simulation_context *sc;
	struct integration_settings settings;
	const char **names;

	parse_args(argc, argv);

	if (!(sc = simulation_context_create_from_sbml_file(model_filename)))
		goto bailout;

	if ((names = simulation_get_value_names(sc)))
	{
		unsigned int i;
		printf("Time");
		for (i=0;names[i];i++)
		{
			printf("\t");
			printf(names[i]);
		}
		printf("\n");
	}

	integration_settings_init(&settings);
	settings.sample_func = sample;
	settings.sample_str_func = sample_strings;
	settings.absolute_error = error;
	settings.relative_error = error;
	settings.time = maxtime;
	settings.steps = sample_steps;
	settings.force_interpreted = force_interpreted;
	settings.stochastic = stochastic;
	settings.stiff = stiff;

	simulation_integrate(sc,&settings);

	simulation_context_free(sc);

	return EXIT_SUCCESS;

bailout:
	return EXIT_FAILURE;
}


