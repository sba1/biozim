#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <config.h>

#define HAVE_LIBGMP

#ifdef HAVE_LIBGMP
#include <gmp.h>
#endif

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

/** @brief The seed was supplied by the user */
static unsigned int user_seed;

/** @brief The seed to be used */
static unsigned int seed;

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

struct variable_assignment *va_first;
struct variable_assignment *va_last;

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


/**
 * @brief indicates whether the mean should be calculated. Only suitable when runs > 1.
 */
static int take_mean;

/**
 * @brief identifies the current run.
 */
static int current_run;



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
			"\t-h, --help                show this help and quit.\n"
			"\t    --error               specifies the error (absolute and relative) given in double.\n"
			"\t    --force-interpreted   forces the interpreted calculation of the rhs.\n"
			"\t    --list-values         list all values.\n"
			"\t    --maxtime             specifies the end time (defaults to 1).\n"
			"\t    --mean                specifies whether the mean should be calculated when runs > 1."
			"\t    --output-time         outputs the current time (stderr).\n"
			"\t    --plot [sp1,...,spn]  plots the results using gnuplot. Optionally, you can specify\n"
			"\t                          the species to be plotted.\n"
			"\t    --runs                specifies the runs to be performed when in stochastic mode.\n"
			"\t    --sample-steps        the number of sample steps (defaults to 5000).\n"
			"\t    --seed                specifies the seed to be used for stochastic simulation.\n"
			"\t    --set-var \"name=dbl\"  set the given variable to double value.\n"
			"\t    --stiff               use solver for stiff ODEs.\n"
			"\t    --stochastic          apply stochastic simulation.\n"
			"\t    --verbose             verbose output.\n",
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
		} else if (!strcmp(argv[i],"--seed"))
		{
			char *nr_arg;

			if (argv[i][6]=='=')
				nr_arg = &argv[i][7];
			else
			{
				nr_arg = argv[++i];
				if (i>=argc)
				{
					fprintf(stderr,"The --seed option needs an integer argument.\n");
					exit(-1);
				}
			}

			seed = strtol(nr_arg,NULL,10);
			user_seed = 1;
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
		} else if (!strcmp(argv[i],"--mean"))
		{
			take_mean = 1;
		} else if (!strcmp(argv[i],"--plot"))
		{
			char *nr_arg;

			plot = 1;

			if (argv[i][6]=='=')
				nr_arg = &argv[i][7];
			else
			{
				if (i + 1 >= argc)
					continue;

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
		} else if (!strcmp(argv[i],"--set-var"))
		{
			struct variable_assignment *sva;
			double value;
			char *nr_arg;
			char *nend;
			char *name;

			if (argv[i][8]=='=')
				nr_arg = &argv[i][9];
			else
			{
				nr_arg = argv[i+1];
				i++;
				if (i>=argc)
				{
					fprintf(stderr,"The --set-var option needs an assignment.\n");
					exit(-1);
				}
			}

			if (!(nend = strchr(nr_arg,'=')))
			{
				fprintf(stderr,"The --set-var option needs an assignment of form \"name=dbl\".\n");
				exit(-1);
			}

			if (!(name = (char*)malloc(nend-nr_arg)))
			{
				fprintf(stderr,"Not enough memory!");
				exit(-1);
			}

			strncpy(name,nr_arg,nend-nr_arg);
			value = strtod(nend+1, NULL);

			if (!(sva = (struct variable_assignment*)malloc(sizeof(*sva))))
			{
				fprintf(stderr,"Not enough memory!");
				exit(-1);
			}

			sva->value_name = name;
			sva->value = value;

			if (!va_first) va_first = sva;
			else va_last->next = sva;
			va_last = sva;
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
 * A single recorded sample
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

static int track_sample(double time, int num_values, double *values)
{
	struct sample *current;
	int i;
	
	if (!(current = (struct sample*)malloc(sizeof(*current))))
		return 0;
	memset(current,0,sizeof(*current));
	
	/* Prepare sample */
	if (!(current->values = (double*)malloc(num_values*sizeof(double))))
		goto bailout;

	current->num_values = num_values;
	for (i=0;i<num_values;i++)
		current->values[i] = values[i];
	current->time = time;

	/* Enqueue sample */
	if (!samples_first) samples_first = current;
	if (samples_last)
		samples_last->next = current;
	samples_last = current;
	samples_total++;

	return 1;
bailout:
	free(current->values);
	free(current);
	return 0;
}



/**
 * Our sampling function.
 *
 * @param time
 * @param num_values
 * @param values
 * @return
 */
int sample(double time, int num_values, double *values, void *user_data)
{
	int i;

	if (output_time)
		fprintf(stderr,"%g\n",time);

	/* Console output */
	printf("%g",time);
	if (runs > 1)
		printf("\t%d",current_run);
	for (i=0;i<num_values;i++)
		printf("\t%.12g",values[i]);

	if (plot)
		track_sample(time,num_values,values);

leave:
	return 1;
}

#ifdef HAVE_LIBGMP

static double *stat_time;
static mpz_t *stat_sum;
static int stat_columns;
static int stat_rows;

/** @brief holds the current sample row. Needs to be reset for a next run */
static int stat_row;

int stat_sample(double time, int num_values, double *values, void *user_data)
{
	if (!stat_sum)
	{
		unsigned int total = (sample_steps+1)*num_values;

		stat_columns = num_values;

		if (!(stat_sum = (mpz_t*)malloc(total * sizeof(mpz_t))))
			return 0;

		for (unsigned int i=0;i<total;i++)
			mpz_init(stat_sum[i]);
	}
	
	if (!stat_time)
	{
		if (!(stat_time = (double*)malloc(sizeof(double)*(sample_steps+1))))
			return 0;

		for (int i=0;i<num_values;i++)
			stat_time[i] = 0.0;
	}


	stat_time[stat_row] = time;
	mpz_t *start = &stat_sum[stat_row * num_values];
	
	for (int i = 0; i < num_values;i++)
		mpz_add_ui(start[i],start[i],values[i]);

	stat_row++;
	if (stat_row > stat_rows) stat_rows = stat_row;
}
#endif

/**
 * Our sampling function for string. This is called always
 * in sync with the above function.
 *
 * @param time
 * @param num_strings
 * @param values
 * @return
 */
int sample_strings(double time, int num_strings, char **values, void *user_data)
{
	int i;

	for (i=0;i<num_strings;i++)
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

	if (!(sc = simulation_context_create_from_sbml_file(model_filename,va_first)))
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
		if (runs > 1)
			printf("\tRun");
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

	if (take_mean)
	{
		settings.sample_func = stat_sample;
		settings.sample_str_func = NULL;
	}

	/* Do the seed stuff */
	if (user_seed)
	{
		settings.seed = seed;
	} else
	{
		unsigned int seed;

		FILE *urand = fopen("/dev/urandom","rb");
		if (urand)
		{
			fread(&seed,1,sizeof(seed),urand);
			fclose(urand);
		} else seed = time(NULL) + clock();

		settings.seed = seed;
	}

	if (runs > 1 && !stochastic)
	{
		fprintf(stderr,"The --runs parameter has been ignored as we don't run in stochastic mode.\n");
		runs = 1;
	}

	for (run=0;run<runs;run++)
	{
		stat_row = 0;
		simulation_context_reset(sc);
		current_run = run;
		simulation_integrate(sc,&settings);
	}

#ifdef HAVE_LIBGMP
	if (take_mean && stat_sum)
	{
		double *values = (double*)malloc(sizeof(double)*stat_columns);
		if (values)
		{
			mpz_t temp;
			mpz_init(temp);

			for (int i=0;i<stat_rows;i++)
			{
				for (int j=0;j<stat_columns;j++)
				{
					unsigned long int r = mpz_fdiv_q_ui(temp,stat_sum[i*stat_columns + j],runs);
					values[j] = mpz_get_ui(temp) + r / (double)runs;
				}
				
				sample(stat_time[i],stat_columns,values,NULL);
				sample_strings(stat_time[i],0,NULL,NULL);
			}
			free(values);
		}
	}
#endif
	

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
			int j;

			for (j=0;plot_species[j];j++)
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
			int idx;

			for (idx=0;idx<samples_first->num_values;idx++)
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
