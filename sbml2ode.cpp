#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Own headers */
#include "simulation.h"

/********************************************************/

static const char *model_filename;
static int force_interpreted;
static int stiff;
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
			"\t    --verbose            verbose output.\n"
			"\t    --force-interpreted  forces the interpreted calculation of the rhs.\n"
			"\t    --stiff              use solver for stiff ODEs.\n",
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
		} else if (!strcmp(argv[i],"--stiff"))
		{
			stiff = 1;
		} else
		{
			filename_given = 1;
			if (strcmp(argv[i],"-"))
				model_filename = argv[i];
		}
	}
	
	if (!filename_given)
	{
		fprintf(stderr,"No filename has been specifed!\n");
		usage(argv[0]);
		exit(-1);		
	}
}

/**********************************************************
 Our sampling function.
***********************************************************/
int sample(double time, int num_values, double *values)
{
	printf("%g",time);
	
	for (int i=0;i<num_values;i++)
	{
		printf("\t%g",values[i]);
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
	char **names;

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
	settings.absolute_error = 1e-15;
	settings.relative_error = 1e-15;
	settings.time = 2000;
	settings.steps = 2000;
	settings.force_interpreted = force_interpreted;
	settings.stiff = stiff;

	simulation_integrate(sc,&settings);

	simulation_context_free(sc);
	
	return EXIT_SUCCESS;
	
bailout:
	return EXIT_FAILURE;
}


