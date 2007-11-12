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
	
	/* The node of ast describing the right part of the ODE (if any) */
	ASTNode *node;
};

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
};


#define environment_set_value_double(v,d) ((struct value*)v)->value = d

void environment_init(struct environment *env, struct environment *parent);
struct value *environment_get_value(const struct environment *env, const char *name);
double environment_get_value_by_name(const struct environment *env, const char *name);
void environment_query_all(const struct environment *env, int (*callback)(struct value *));
int environment_is_value_defined(const struct environment *env, const char *name);
void *environment_get_value_handle(const struct environment *env, const char *name);
struct value *environment_add_value(struct environment *env, char *name);

#endif
