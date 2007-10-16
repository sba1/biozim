#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Own headers */
#include "simulation.h"

/********************************************************/

static const char *model_filename;
int verbose; /* used by other modules */

/*************************************************
 Displays the program's usage and exits
*************************************************/
static void usage(char *name)
{
	fprintf(stderr, "usage: %s [OPTIONS] SBML-File ...\n"
			"Loads the given SBML-File and performs a simulation run.\n"
			"Specify '-' to read from stdin.\n"
			"\t-h, --help      show this help and quit.\n",
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
		if (!strcmp(argv[i],"-h") && !strcmp(argv[i],"--help"))
		{
			usage(argv[0]);
			exit(-1);
		} else if (!strcmp(argv[i],"--verbose"))
		{
			verbose = 1;
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
 Main Entry
***********************************************************/
int main(int argc, char **argv)
{
	struct simulation_context *sc;
	struct integration_settings settings;

	parse_args(argc, argv);

	if (!(sc = simulation_context_create_from_sbml_file(model_filename)))
		goto bailout;

	integration_settings_init(&settings);
	
	simulation_integrate(sc,&settings);

	simulation_context_free(sc);
	

	// Now create settings object for integration of model

/*	ODESettings *settings = new ODESettings();
	cout << *settings;
	// Update some settings
	double time = 1000.0;
	int printstep = 100;
	settings->setTimePointSeries(time, printstep);
	settings->setAbsoluteError(1e-9);
	settings->setRelativeError(1e-4);
	settings->setMaximumSteps(1000);
	cout << *settings;
*/

	return EXIT_SUCCESS;
	
bailout:
	return EXIT_FAILURE;
}


