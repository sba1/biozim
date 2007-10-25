#ifndef SIMULATION_H_
#define SIMULATION_H_

/* Opaque */
struct simulation_context;

/* Data structure for settings */
struct integration_settings
{
	/* Specifies the end time starting from 0.0 */
	double time;
	
	/* Number of sampling steps within the time */
	unsigned int steps;

	double absolute_error;
	double relative_error;

	/* The sampling function. It is called for every sampled time point */
	int (*sample_func)(double time, int num_values, double *values);

	/* If set, the interpeted evaluation of the right-hand side is performed */
	int force_interpreted;
	
	/* Use integerator for stiff systems */
	int stiff;
};

void integration_settings_init(struct integration_settings *settings);

struct simulation_context *simulation_context_create_from_sbml_file(const char *filename);
char ** simulation_get_value_names(struct simulation_context *sc);
void simulation_integrate(struct simulation_context *sc, struct integration_settings *settings);
void simulation_context_free(struct simulation_context *sc);

#endif /*SIMULATION_H_*/
