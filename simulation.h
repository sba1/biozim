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

	/* As above, but for string variables. The calling of both functions is synced */
	int (*sample_str_func)(double time, int num_strings, char **values);

	/* If set, the interpeted evaluation of the right-hand side is performed */
	int force_interpreted;

	/* Use integerator for stiff systems */
	int stiff;

	/* Use stochastic simulator */
	int stochastic;
};

void integration_settings_init(struct integration_settings *settings);

struct simulation_context *simulation_context_create_from_sbml_file(const char *filename);
void simulation_context_reset(struct simulation_context *sc);
const char ** simulation_get_value_names(struct simulation_context *sc);
int simulation_context_get_value_index(struct simulation_context *sc, const char *name);
void simulation_integrate(struct simulation_context *sc, struct integration_settings *settings);
void simulation_context_free(struct simulation_context *sc);

#endif /*SIMULATION_H_*/
