#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Own headers */
#include "simulation.h"
#include "gnuplot_i.h"

/********************************************************/

/** @brief the filename of the sbml file to be processed */
static const char *model_filename;

/** @brief don't use JIT compiling technique to speed up the simulation */
static int no_jit;

/**
 * @brief Use stiff the ODE solver
 *
 * @note Only applicable for deterministic solver
 */
static int stiff;

/** @brief Perform stochastic simulation */
static int stochastic;

/** @brief Print time to stderr */
static int output_time;

/** @brief Specifies the time point at which the simulation is stopped */
static double maxtime;

/** @brief List values */
static int list_values;

/** @brief Specifies the error value to be used in deterministic simulation */
static double error;

/** @brief Number of total samples to be taken */
static int sample_steps;

/** @brief Species if the result should be plotted */
static int plot;

/**
 * @brief Species the species to be considered in the plot
 *
 * Can be NULL to indicate that all species needs to be plotted
 */
static char **plot_species;

/**
 * @brief Specifies the number of runs to be performed
 *
 * @note Values greater than 1 make sense only when simulation
 * runs in stochastic mode
 */
static int runs = 1;

/** @brief Verbose output to stderr */
int verbose; /* used by other modules */


/********************************************************/

/**
 * Displays the program's usage and exits.
 *
 * @param name
 */
static void usage(char *name)
{
	fprintf(stderr, "usage: %s [OPTIONS] SBML-File ...\n"
			"Loads the given SBML-File and performs a simulation run.\n"
			"Specify '-' to read from stdin.\n"
			"\t-h, --help               show this help and quit.\n"
			"\t    --error              specifies the error (absolute and relative) given in double.\n"
			"\t    --force-interpreted  forces the interpreted calculation of the rhs.\n"
			"\t    --list-values        list all values.\n"
			"\t    --maxtime            specifies the end time (defaults to 1).\n"
			"\t    --output-time        outputs the current time (stderr).\n"
			"\t    --plot [sp1,...,spn] plots the results using gnuplot. Optionally, you can specify\n"
			"\t                         the species to be plotted.\n"
			"\t    --runs               specifies the runs to be performed when in stochastic mode.\n"
			"\t    --sample-steps       the number of sample steps (defaults to 5000).\n"
			"\t    --stiff              use solver for stiff ODEs.\n"
			"\t    --stochastic         apply stochastic simulation.\n"
			"\t    --verbose            verbose output.\n",
			name);

	exit(1);
}


/**
 * Splits the given character sequence.
 *
 * @param str
 * @param c
 * @return
 */
static char **strsplit(char *str, char c)
{
	int count, pos;
	char **array;
	char *dup;

	if (!(dup = strdup(str)))
		return NULL;

	pos = count = 0;
	while (dup[pos])
	{
		if (dup[pos] == c)
			count++;
		pos++;
	}

	if (!(array = (char**)malloc(sizeof(char*) * (count + 2))))
		return NULL;

	array[0] = dup;
	array[count+1] = NULL;

	pos = count = 0;
	while (dup[pos])
	{
		if (dup[pos] == c)
		{
			array[++count] = &dup[pos+1];
			dup[pos] = 0;
		}
		pos++;
	}

	return array;
}

/**
 * Returns the index of the given string (or -1 if string
 * couldn't be found) in the given NULL-terminated array.
 *
 * @param array
 * @param str
 * @return
 */
static int strindex(const char **array, const char *str)
{
	int i = 0;

	while (array[i])
	{
		if (!strcmp(array[i],str))
			return i;
		i++;
	}
	return -1;
}

/**
 * Parses the command line.
 *
 * TODO: When increasing i, check for out of bounds error.
 *
 * @param argc
 * @param argv
 */
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
		} else if (!strcmp(argv[i],"--list-values"))
		{
			list_values = 1;
		} else if (!strcmp(argv[i],"--force-interpreted"))
		{
			no_jit = 1;
		} else if (!strcmp(argv[i],"--runs"))
		{
			char *nr_arg;

			if (argv[i][7]=='=')
				nr_arg = &argv[i][8];
			else
			{
				nr_arg = argv[++i];
				if (i>=argc)
				{
					fprintf(stderr,"The --sample-steps option needs an integer argument.\n");
					exit(-1);
				}
			}

			runs = strtod(nr_arg, NULL);
			if (runs < 0) runs = 0;
		} else if (!strcmp(argv[i],"--plot"))
		{
			char *nr_arg;

			plot = 1;

			if (argv[i][6]=='=')
				nr_arg = &argv[i][7];
			else
			{
				nr_arg = argv[i+1];
				if (*nr_arg == '-')
					continue;
				i++;
			}
			if (!(plot_species = strsplit(nr_arg,',')))
				fprintf(stderr,"Could not determine the name of species to be plotted!");
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

			if (argv[i][14]=='=')
				nr_arg = &argv[i][15];
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

			if (argv[i][7]=='=')
				nr_arg = &argv[i][8];
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
			if (strcmp(argv[i],"-"))
			{
				model_filename = argv[i];
				filename_given = 1;
			} else
			{
				fprintf(stderr,"Unknown parameter \"%s\"\n",argv[i]);
				usage(argv[0]);
				exit(-1);
			}
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

/**
 * A single sample
 */
struct sample
{
	struct sample *next;

	double time;
	int num_values;
	double *values;
	int num_strings;
	char **strings;
};

static struct sample *samples_first;
static struct sample *samples_last;
static int samples_total;

/**
 * Our sampling function.
 *
 * @param time
 * @param num_values
 * @param values
 * @return
 */
int sample(double time, int num_values, double *values)
{
	if (output_time)
		fprintf(stderr,"%g\n",time);

	/* Console output */
	printf("%g",time);
	for (int i=0;i<num_values;i++)
		printf("\t%.12g",values[i]);

	if (plot)
	{
		struct sample *current = (struct sample*)malloc(sizeof(*current));
		if (current)
		{
			/* Prepare sample */
			memset(current,0,sizeof(*current));
			if (!(current->values = (double*)malloc(num_values*sizeof(double))))
			{
				free(current);
				goto leave;
			}

			current->num_values = num_values;
			for (int i=0;i<num_values;i++)
				current->values[i] = values[i];
			current->time = time;

			/* Enqueue sample */
			if (!samples_first) samples_first = current;
			if (samples_last)
				samples_last->next = current;
			samples_last = current;
			samples_total++;
		}
	}

leave:
	return 1;
}

/**
 * Our sampling function for string. This is called always
 * in sync with the above function.
 *
 * @param time
 * @param num_strings
 * @param values
 * @return
 */
int sample_strings(double time, int num_strings, char **values)
{
	for (int i=0;i<num_strings;i++)
	{
		printf("\t\"%s\"",values[i]);
	}
	printf("\n");

	return 1;
}

/**
 * Plots the given value.
 *
 * @param ctrl
 * @param which
 * @return
 */
void plot_samples(gnuplot_ctrl *ctrl, int which, const char *name, double *x, double *y)
{
	struct sample *s = samples_first;
	int i = 0;

	while (s)
	{
		x[i] = s->time;
		y[i] = s->values[which];

		i++;
		s = s->next;
	}
	gnuplot_plot_xy(ctrl,x,y,samples_total,name);
}

#include "environment.h"

/**
 * Print value.
 *
 * @param v
 * @return
 */
static int print_value(struct value *v)
{
	printf("%s\t%g\n",v->name,v->value);
	return 1;
}

/**
 * Main entry
 *
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv)
{
	struct simulation_context *sc;
	struct integration_settings settings;
	const char **names;
	int run;

	parse_args(argc, argv);

	if (!(sc = simulation_context_create_from_sbml_file(model_filename)))
		goto bailout;

	names = simulation_get_value_names(sc);

	if (list_values)
	{
		if (names)
			simulation_context_query_values(sc,print_value,NULL);
		simulation_context_free(sc);
		return 0;
	}

	if (names)
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
	settings.force_interpreted = no_jit;
	settings.stochastic = stochastic;
	settings.stiff = stiff;

	if (runs > 1 && !stochastic)
	{
		fprintf(stderr,"The --runs parameter has been ignored as we don't run in stochastic mode.\n");
		runs = 1;
	}

	for (run=0;run<runs;run++)
	{
		simulation_context_reset(sc);
		simulation_integrate(sc,&settings);
	}

	if (plot && samples_first)
	{
		double *x, *y;
		gnuplot_ctrl *ctrl;

		if (!(x = (double*)malloc(samples_total*sizeof(double))))
			goto bailout;
		if (!(y = (double*)malloc(samples_total*sizeof(double))))
			goto bailout;

		if (!(ctrl = gnuplot_init()))
		{
			fprintf(stderr,"Couldn't initialize gnuplot!");
			goto bailout;
		}

		if (plot_species)
		{
			for (int j=0;plot_species[j];j++)
			{
				int idx = strindex(names,plot_species[j]);
				if (idx == -1)
				{
					fprintf(stderr,"Species \"%s\" couldn't not be found!\n",plot_species[j]);
					continue;
				}
				plot_samples(ctrl,idx,names[idx],x,y);
			}
		} else
		{
			for (int idx=0;idx<samples_first->num_values;idx++)
				plot_samples(ctrl,idx,names[idx],x,y);
		}

		fprintf(stderr,"Press Enter to quit!\n");
		while (getchar()!='\n') {}
		gnuplot_close(ctrl);
	}

	simulation_context_free(sc);

	return EXIT_SUCCESS;

bailout:
	return EXIT_FAILURE;
}
