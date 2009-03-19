#ifndef __ENVIRONMENT_H
#define __ENVIRONMENT_H

class ASTNode;

struct value
{
	/* The name of the value */
	char *name;

	/* The value's actual value (used for ODEs) */
	double value;

	/* molecules (used for stochastic simulation) */
	int molecules;

	/* Whether fixed */
	int fixed;

	/* Contains the index in the environment */
	int index;

	/* Variable represents a species */
	int is_species;

	/* Variable is uninitialized */
	int uninitialized;

	int has_only_substance_units;
	struct value *compartment_value;

	/* The node of ast describing the right part of the ODE (if any) */
	ASTNode *node;
};

struct environment_internal;

/* Container for values */
struct environment
{
	/* NULL for the parent environment */
	struct environment *parent;

	/* Array of values */
	struct value **values;

	/* Length of the values array */
	unsigned int num_values;

	/* Allocated length of values */
	unsigned int num_values_allocated;

	/* Internal stuff */
	struct environment_internal *internal;
};

struct environment_snapshot
{
	int num_entries;
	double *dbl_values;
	int *int_values;
};

#define environment_set_value_double(v,d) ((struct value*)v)->value = d

void environment_init(struct environment *env, struct environment *parent);
void environment_deinit(struct environment *env);

void environment_optimize(struct environment *env);
struct environment_snapshot *environment_snapshot(struct environment *env);
void environment_set_to_snapshot(struct environment *env, struct environment_snapshot *snap);
void environment_free_snapshot(struct environment_snapshot *snap);
struct value *environment_get_value(const struct environment *env, const char *name);
double environment_get_value_by_name(const struct environment *env, const char *name);
void environment_query_all(const struct environment *env, int (*callback)(struct value *));
int environment_is_value_defined(const struct environment *env, const char *name);
void *environment_get_value_handle(const struct environment *env, const char *name);
struct value *environment_add_value(struct environment *env, const char *name);

#endif
