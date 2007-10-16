#ifndef SIMULATION_H_
#define SIMULATION_H_

/* Opaque */
struct simulation_context;

struct simulation_context *simulation_context_create_from_sbml_file(const char *filename);
void simulation_integrate(struct simulation_context *sc);
void simulation_context_free(struct simulation_context *sc);

#endif /*SIMULATION_H_*/
