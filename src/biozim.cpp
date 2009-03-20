#include <limits.h>
#include <math.h>
#include <search.h>
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

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

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

/** @brief Specifies if the result should be plotted */
static int plot;

/** @brief Specifies if the distribution should be plotted */
static int plot_distribution;

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
 * @brief indicates whether the standard deviation should be calculated. Only suitable when runs > 1.
 */
static int take_stddev;

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
			"\t    --plot-distribution   plots the result as distribution.\n"
			"\t    --runs                specifies the runs to be performed when in stochastic mode.\n"
			"\t    --sample-steps        the number of sample steps (defaults to 5000).\n"
			"\t    --stddev              specifies whether the standard deviation should be calculated when runs > 1.\n"
			"\t    --seed                specifies the seed to be used for stochastic simulation.\n"
			"\t    --set-var \"name=dbl\"  set the given variable to double value.\n"
			"\t    --stiff               use solver for stiff ODEs.\n"
			"\t    --stochastic          apply stochastic simulation.\n"
			"\t    --verbose             verbose output.\n",
			name);

	exit(1);
}


/**
 * Splits the given character sequence. The returned array is NULL terminated.
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
static int strindex(const char * const *array, const char *str)
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
		} else if (!strcmp(argv[i],"--stddev"))
		{
			take_stddev = 1;
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
		} else if (!strcmp(argv[i],"--plot-distribution"))
		{
			plot_distribution = 1;
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
static int current_slot;
static int *samples_max_values;
static int *samples_min_values;

static void **hist_vector; /* every element contains a binary tree corresponding to a time slot and a variable */
static void ***hist; /* time slots are rows, values are columns */
static int *max_count;

struct hist_item
{
	int value;
	int count;
};

int hist_item_compare(const void *pa, const void *pb)
{
	const struct hist_item *h1 = (struct hist_item*)pa;
	const struct hist_item *h2 = (struct hist_item*)pb;
	
//	printf("%p(%d) %p(%d)\n",h1,h1->value,h2,h2->value);
	return h1->value - h2->value;
}

static int track_sample(double time, int num_values, double *values)
{
	struct sample *current;
	int i;

	/* Allocate memory when called for the first time */
	if (!samples_min_values)
	{
		if ((samples_min_values = (int*)malloc(num_values * sizeof(int))))
		{
			for (i=0;i<num_values;i++)
				samples_min_values[i] = INT_MAX;
		}
	}
	if (!samples_max_values)
	{
		if ((samples_max_values = (int*)malloc(num_values * sizeof(int))))
		{
			for (i=0;i<num_values;i++)
				samples_min_values[i] = 0;
		}
	}

	if (!hist_vector && plot_distribution)
	{
		if ((hist_vector = (void**)malloc((sample_steps+1)*num_values*sizeof(hist_vector[0]))))
		{
			memset(hist_vector,0,(sample_steps+1)*num_values*sizeof(hist_vector[0]));
			if ((hist = (void***)malloc((sample_steps+1)*sizeof(hist[0]))))
			{
				for (i=0;i<sample_steps+1;i++)
					hist[i] = &hist_vector[i*num_values];
			}
		}
	}

	if (!max_count && plot_distribution)
	{
		if ((max_count = (int*)malloc(num_values*sizeof(int))))
			memset(max_count,0,num_values*sizeof(int));
	}

	if (!(current = (struct sample*)malloc(sizeof(*current))))
		return 0;
	memset(current,0,sizeof(*current));
	
	/* Prepare sample */
	if (!(current->values = (double*)malloc(num_values*sizeof(double))))
		goto bailout;

	current->num_values = num_values;
	for (i=0;i<num_values;i++)
		current->values[i] = values[i];

	if (samples_min_values)
		for (i=0;i<num_values;i++)
			samples_min_values[i] = MIN(values[i],samples_min_values[i]);

	if (samples_max_values)
		for (i=0;i<num_values;i++)
			samples_max_values[i] = MAX(values[i],samples_max_values[i]);

	if (hist)
	{
		for (i=0;i<num_values;i++)
		{
			struct hist_item item;
			item.value = values[i];
			item.count = 1;

			struct hist_item **item_found;
			
			if (!(item_found = (struct hist_item**)tfind(&item,&hist[current_slot][i],hist_item_compare)))
			{
				struct hist_item *new_item;
				
				if ((new_item = (struct hist_item*)malloc(sizeof(struct hist_item))))
				{
					new_item->count = 1;
					new_item->value = item.value;
					
					if (!tsearch(new_item,&hist[current_slot][i],hist_item_compare))
						fprintf(stderr,"Unable to add an item\n");
				}
			} else
			{
				(*item_found)->count++;
				
				if (max_count)
				{
					if ((*item_found)->count > max_count[i])
						max_count[i] = (*item_found)->count; 
				}
			}
		}
	}
	current->time = time;

	/* Enqueue sample */
	if (!samples_first) samples_first = current;
	if (samples_last)
		samples_last->next = current;
	samples_last = current;
	samples_total++;
	current_slot++;

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
	if (runs > 1 && (take_mean + take_stddev) != 1)
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
static mpz_t *stat_sum_of_sqr;
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
	
	if (!stat_sum_of_sqr && take_stddev)
	{
		unsigned int total = (sample_steps+1)*num_values;

		if (!(stat_sum_of_sqr = (mpz_t*)malloc(total * sizeof(mpz_t))))
			return 0;

		for (unsigned int i=0;i<total;i++)
			mpz_init(stat_sum_of_sqr[i]);
	}
	
	if (!stat_time)
	{
		if (!(stat_time = (double*)malloc(sizeof(double)*(sample_steps+1))))
			return 0;

		for (int i=0;i<=sample_steps;i++)
			stat_time[i] = 0.0;
	}

	if (stat_time[stat_row] != 0.0 && stat_time[stat_row] != time)
		fprintf(stderr,"***Warning***: Time slices don't match (%g != %g)\n",stat_time[stat_row],time);
	
	stat_time[stat_row] = time;

	mpz_t *start_sum = &stat_sum[stat_row * num_values];
	for (int i = 0; i < num_values;i++)
		mpz_add_ui(start_sum[i],start_sum[i],values[i]);
	
	if (take_stddev)
	{
		mpz_t temp;
		mpz_init(temp);
	
		mpz_t *start_sum_of_sqr = &stat_sum_of_sqr[stat_row * num_values];
		for (int i = 0; i < num_values;i++)
		{
			mpz_set_ui(temp,values[i]);
			mpz_mul_ui(temp,temp,values[i]);
			mpz_add(start_sum_of_sqr[i],start_sum_of_sqr[i],temp);
		}
	}
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

static int hist_item_count = 0;

static void hist_item_count_action(const void *nodep, const VISIT which, const int depth)
{
	struct hist_item *data;

	switch (which)
	{
	case	preorder:
			break;
	case	endorder:
			break;
	case	postorder:
	case	leaf:
			hist_item_count++;
			break;
	}
}

static int hist_item_current;
static struct hist_item **hist_item_array;

static void hist_item_linerize_action(const void *nodep, const VISIT which, const int depth)
{
	struct hist_item *data;

	switch (which)
	{
	case	preorder:
			break;
	case	endorder:
			break;
	case	postorder:
	case	leaf:
			data = *(struct hist_item**)nodep;
			hist_item_array[hist_item_current++] = data;
			break;
	}
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

	parse_args(argc, argv);

	if (runs < 2 && take_mean)
	{
		fprintf(stderr,"The mean has been requested but the number of requested runs (--run parameter) wasn't higher than 1!\n");
		goto bailout;
	}

	if (runs < 2 && take_stddev)
	{
		fprintf(stderr,"The standard deviation has been requested but the number of requested runs (--run parameter) wasn't higher than 1!\n");
		goto bailout;
	}

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
		
		if (runs > 1 && (take_mean + take_stddev) != 1)
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

	if (take_mean || take_stddev)
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

	for (current_run=0;current_run<runs;current_run++)
	{
		stat_row = 0;
		current_slot = 0;
		simulation_context_reset(sc);
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

			current_run = 0;

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


	if (take_stddev && stat_sum && stat_sum_of_sqr)
	{
		double *values = (double*)malloc(sizeof(double)*stat_columns);
		if (values)
		{
			mpz_t sqr_sum;
			mpz_t temp;
			mpq_t var;
 			mpq_t divisor; /* will contain n * (n-1) */

			mpz_init(sqr_sum);
			mpz_init(temp);
			mpq_init(var);
			mpq_init(divisor);
			
			mpz_set_ui(temp,runs);
			mpz_mul_ui(temp, temp, runs - 1);
			mpq_set_z(divisor,temp);

			current_run = 1;

			for (int i=0;i<stat_rows;i++)
			{
				for (int j=0;j<stat_columns;j++)
				{
					mpz_set(sqr_sum, stat_sum[i*stat_columns + j]);          /* sqr_sum = sum */
					mpz_mul(sqr_sum, sqr_sum, stat_sum[i*stat_columns + j]); /* sqr_sum = sum^2 */ 

					mpz_set(temp, stat_sum_of_sqr[i*stat_columns + j]);		 /* temp = sum_of_sqr */
					mpz_mul_ui(temp, temp, runs);							 /* temp = n * sum_of_sqr */
					mpz_sub(temp,temp,sqr_sum);								 /* temp = n * sum_of_sqr - sum^2 */

					mpq_set_z(var,temp);									 /* var = n * sum_of_sqr - sum^2 */ 
					mpq_div(var,var,divisor);								 /* var = (n * sum_of_sqr - sum^2) / (n*(n-1)) */
					
					values[j] = sqrt(mpq_get_d(var));
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

		if (plot_distribution)
		{
			for (int i=0;i<samples_first->num_values;i++)
			{
				if (plot_species && strindex(plot_species,names[i]) == -1) continue;

				int obj_index = 2;

				gnuplot_cmd(ctrl,"set object 1 rect at 0,0 size %d,%d fc rgbcolor \"#000000\" lw 0\n",current_slot*2,samples_max_values[i]*2);

				/* current slot holds the last time slot */
				for (int j=0;j<current_slot;j++)
				{
					void *root = hist[j][i];
					
					hist_item_count = 0;
					twalk(root,hist_item_count_action);

					free(hist_item_array);
					if ((hist_item_array = (struct hist_item**)malloc(hist_item_count * sizeof(struct hist_item*))))
					{
						hist_item_current = 0;
						twalk(root,hist_item_linerize_action);

						for (int k=0;k<hist_item_count;k++)
						{
							struct hist_item *hi = hist_item_array[k];
							double x = j;
							double y = hi->value;

							int r = 255 * (log(hi->count) / log(max_count[i]));
							
							int g = 0;
							int b = 0;
							
//							fprintf(stderr,"set object %d rect at %g,%g size 1.1,1.1 fc rgbcolor \"#%02x%02x%02x\" lw 0\n",obj_index,x+0.5,y+0.5,r,g,b);
							gnuplot_cmd(ctrl,"set object %d rect at %g,%g size 1.1,1.1 fc rgbcolor \"#%02x%02x%02x\" lw 0\n",obj_index,x+0.5,y+0.5,r,g,b);
							obj_index++;
						}
					}
				}
				
				gnuplot_cmd(ctrl,"plot [0:%d] [0:%d] 0",current_slot,samples_max_values[i]);
			}
		} else
		{
			if (plot_species)
			{
				int j;
	
				for (j=0;plot_species[j];j++)
				{
					int idx = strindex(names,plot_species[j]);
					if (idx == -1)
					{
						fprintf(stderr,"Species \"%s\" couldn't be found!\n",plot_species[j]);
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
