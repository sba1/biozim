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
};

void integration_settings_init(struct integration_settings *settings);

struct simulation_context *simulation_context_create_from_sbml_file(const char *filename);
void simulation_integrate(struct simulation_context *sc, struct integration_settings *settings);
void simulation_context_free(struct simulation_context *sc);

#endif /*SIMULATION_H_*/