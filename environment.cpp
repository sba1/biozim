#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "environment.h"

/*****************************************************
 Initializes the environment
******************************************************/
void environment_init(struct environment *env, struct environment *parent)
{
	memset(env,0,sizeof(*env));
	env->parent = parent;
}

/*****************************************************
 Returns the reference to the value or NULL, if
 doesn't exist.
******************************************************/
struct value *environment_get_value(const struct environment *env, const char *name)
{
	unsigned int i;

	while (env)
	{
		for (i=0;i<env->num_values;i++)
		{
			if (!strcmp(name,env->values[i]->name))
				return env->values[i];
		}
		env = env->parent;
	}
	return NULL;
}

/*****************************************************
 Returns the value of a varibale by its name.
******************************************************/
double environment_get_value_by_name(const struct environment *env, const char *name)
{
	struct value *v = environment_get_value(env, name);
	if (v) return v->value;
	fprintf(stderr,"***Warning***: Symbol \"%s\" not found! Defaulting to 0.0\n",name);
	return 0.0;
}

/*****************************************************
 Querys the names
******************************************************/
void environment_query_all(const struct environment *env, int (*callback)(struct value *))
{
	unsigned int i;

	for (i=0;i<env->num_values;i++)
	{
		if (!(callback(env->values[i])))
			return;
	}
}

/*****************************************************
 Returns whether the environment variable is defined
 not.
******************************************************/
int environment_is_value_defined(const struct environment *env, const char *name)
{
	return !!environment_get_value(env, name);
}

/*****************************************************
 Returns a handle to the value of the given
 name.
******************************************************/
void *environment_get_value_handle(const struct environment *env, const char *name)
{
	return environment_get_value_handle(env,name);
}

/*****************************************************
 Add a new value to the environment. Returns NULL
 on an error, i.e., if no memory was available.
******************************************************/
struct value *environment_add_value(struct environment *env, const char *name)
{
	struct value *v;

	if (environment_is_value_defined(env,name))
	{
		fprintf(stderr,"Value \"%s\" is already defined!\n",name);
		return NULL;
	}

	if (!(v = (struct value*)malloc(sizeof(struct value))))
	{
		fprintf(stderr,"Not enough memory!\n");
		return NULL;
	}

	if (env->num_values >= env->num_values_allocated)
	{
		struct value **va;

		if (!(va = (struct value**)realloc(env->values,(env->num_values_allocated+10)*2*sizeof(struct value*))))
		{
			fprintf(stderr,"Not enough memory!\n");
			return NULL;
		}
		env->num_values_allocated = (env->num_values_allocated + 10) * 2;
		env->values = va;
	}

	memset(v,0,sizeof(struct value));
	v->name = name;
	v->index = env->num_values;
	env->values[env->num_values++] = v;
	return v;
}

/**
 * Snapshots the current environment.
 *
 * @param env
 * @return
 */
struct environment_snapshot *environment_snapshot(struct environment *env)
{
	struct environment_snapshot *snap;

	unsigned int i;

	if (!(snap = (struct environment_snapshot*)malloc(sizeof(*snap))))
		return NULL;
	if (!(snap->dbl_values = (double*) malloc(env->num_values * sizeof(double))))
	{
		free(snap);
		return NULL;
	}
	if (!(snap->int_values = (int*)malloc(env->num_values * sizeof(int))))
	{
		free(snap->dbl_values);
		free(snap);
		return NULL;
	}

	snap->num_entries = env->num_values;
	for (i=0;i<env->num_values;i++)
	{
		snap->dbl_values[i] = env->values[i]->value;
		snap->int_values[i] = env->values[i]->molecules;
	}
	return snap;
}


/**
 * Sets the values to the content of the snapshot.
 *
 * @param env
 * @param snap
 */
void environment_set_to_snapshot(struct environment *env, struct environment_snapshot *snap)
{
	int i;

	for (i=0;i<snap->num_entries;i++)
	{
		env->values[i]->value = snap->dbl_values[i];
		env->values[i]->molecules = snap->int_values[i];
	}
}

/**
 * Frees all resources associated with the given snapshot.
 *
 * @param snap
 */
void environment_free_snapshot(struct environment_snapshot *snap)
{
	free(snap->int_values);
	free(snap->dbl_values);
	free(snap);
}
