//#define HAVE_CMPH_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef HAVE_CMPH_H
#include <cmph.h>
#endif

#include "environment.h"

struct environment_internal
{
	int dummy;
#ifdef HAVE_CMPH_H
	int names_length;
	char **names;
	cmph_io_adapter_t *cmph_source;
	cmph_t *cmph_hash;
#endif
};

/**
 * Initializes the environment.
 *
 * @param env the environment to be initialized.
 * @param parent the parent environment or NULL.
 */
void environment_init(struct environment *env, struct environment *parent)
{
	memset(env,0,sizeof(*env));
	env->parent = parent;
}

/**
 * Cleanup resources allocated by environment_optimize().
 * 
 * @param env
 */
static void environment_clean_optimization(struct environment_internal *internal)
{
	if (!internal) return;
#ifdef HAVE_CMPH_H
	cmph_io_vector_adapter_destroy(internal->cmph_source);
	if (internal->names)
	{
		int i;
		
		for (i=0;i<internal->names_length;i++)
			free(internal->names[i]);
		free(internal->names);
	}
#endif
	free(internal);
}

/**
 * Optimize the environment for quick accesses. You should
 * call this function when you are read with adding variables.
 *
 * @param env the environment to be optimized.
 */
void environment_optimize(struct environment *env)
{
#ifdef HAVE_CMPH_H
	unsigned int i;
	struct environment_internal *internal;
	cmph_config_t *cmph_config = NULL;

	environment_clean_optimization(env->internal);
	env->internal = NULL;

	/* Allocate new structure */
	if (!(internal = (struct environment_internal *)malloc(sizeof(struct environment_internal))))
		goto bailout;
	memset(internal,0,sizeof(struct environment_internal));

	/* Get names */
	if (!(internal->names = (char**)malloc(env->num_values*sizeof(char*))))
		return;
	for (i=0;i<env->num_values;i++)
	{
		if (!(internal->names[internal->names_length] = strdup(env->values[i]->name)))
			goto bailout;
		internal->names_length++;
	}

	/* Generate perfect hashing */
	if (!(internal->cmph_source = cmph_io_vector_adapter(internal->names,env->num_values)))
		goto bailout;

	if (!(cmph_config = cmph_config_new(internal->cmph_source)))
		goto bailout;

	if (!(internal->cmph_hash = cmph_new(cmph_config)))
		goto bailout;

	cmph_config_destroy(cmph_config);
	env->internal = internal;
	return;

bailout:
	fprintf(stderr,"***Warning***: Couldn't create hash table!");
	environment_clean_optimization(internal);
#endif
}

/**
 * Returns the reference to the value or NULL, if
 * doesn't exist.
 *
 * @param env the environment in which name should be looked up.
 * @param name
 * @return
 */
struct value *environment_get_value(const struct environment *env, const char *name)
{
	unsigned int i;

#ifdef HAVE_CMPH_H	
	struct environment_internal *internal;

	if ((internal = env->internal))
	{
		unsigned int id;

		id = cmph_search(internal->cmph_hash,name,strlen(name));

		if (id < env->num_values && !strcmp(env->values[id]->name,name))
			return env->values[id];
	}
#endif

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
 Returns the value of a variable by its name.
******************************************************/
double environment_get_value_by_name(const struct environment *env, const char *name)
{
	struct value *v = environment_get_value(env, name);
	if (v)
	{
		if (v->uninitialized)
		{
			fprintf(stderr,"***Warning***: Accessing an uninitialized variable \"%s\"!\n",name);
			return 0.0;
		}
		return v->value;
	}
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
		fprintf(stderr,"***Warning***: Value \"%s\" is already defined!\n",name);
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
