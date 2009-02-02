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

	/* If set, the interpreted evaluation of the right-hand side is performed */
	int force_interpreted;

	/* Use integrator for stiff systems */
	int stiff;

	/* Use stochastic simulator */
	int stochastic;
};

/**
 * Thats a single variable assignment.
 */
struct variable_assignment
{
	char *value_name;
	double value;

	struct variable_assignment *next;
};

void integration_settings_init(struct integration_settings *settings);

struct simulation_context *simulation_context_create_from_sbml_file(const char *filename, struct variable_assignment *assignment);
void simulation_context_reset(struct simulation_context *sc);
void simulation_context_set_value(const char *name, double value);
void simulation_context_query_values(struct simulation_context *sc, int (*callback)(struct value *), void *userdata);
const char ** simulation_get_value_names(struct simulation_context *sc);
void simulation_integrate(struct simulation_context *sc, struct integration_settings *settings);
void simulation_context_free(struct simulation_context *sc);

#endif /*SIMULATION_H_*/
