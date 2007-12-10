/* simulation.cpp */

#include <dlfcn.h>
#include <stdint.h>

#include <sbml/SBMLTypes.h>

/* Headers of sundials */
#include <nvector/nvector_serial.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

/* Own headers */
#include "sbml/SBMLParser.h"
#include "ode/ODESettings.h"

#include "environment.h"
#include "simulation.h"

/***********************************************/

struct simulation_context
{
	/* The global environment, here species and global parameters are added */
	struct environment global_env;

	/* Convenience-array of value names, NULL terminated */
	char **names;

	/* Contains indices to values which have an ODE attached */
	struct value **unfixed;

	/* Length of the unfixed array */
	unsigned int num_unfixed;

	/* Containes indiced to values which don't have an ODE attached */
	struct value **fixed;

	/* Length of the fixed array */
	unsigned int num_fixed;
	
	/* All events of this model */
	struct event **events;
	
	/* The size of the events array */
	unsigned int num_events;

	/* Indicates the event's trigger state */
	int *events_active;

	/* Number of reactions */
	unsigned int num_reactions;
	
	/* The reactions */
	struct reaction *reactions;

	/* Dynamic loading support */
	void *dlhandle;
	int (*dlrhs)(double t, double *y, double *ydot, void *f_data);
	int (*dlcheck_trigger)(double t, double *y, int *events_active);

	/* Preallocated space for the dlrhs function */
	double *dly;
	double *dlydot;
};

extern int verbose;

/***********************************************/

/* An event assignment */
struct assignment
{
	struct value *value;
	ASTNode *math;
};

/* An event */
struct event
{
	ASTNode *trigger;

	unsigned int num_assignments;
	struct assignment assignments[0];
};

/* Reference */
struct reference
{
	struct value *value;
	ASTNode *stoich;
};

/* A reaction */
struct reaction
{
	/* For the stochastic simulation */
	double h,c,a;
	ASTNode *formula;

	unsigned int num_reactants;
	struct reference *reactants;

	unsigned int num_products;
	struct reference *products;
};

/***********************************************/


/*****************************************************
 Add a new parameter to the given value.
******************************************************/
static int environment_add_parameter(struct environment *env, Parameter *p)
{
	char *name;
	struct value *v;

	if (!(name = strdup(p->getId().c_str())))
	{
		fprintf(stderr,"Not enough memory!\n");
		return 0;
	}

	if ((v = environment_add_value(env,name)))
	{
		v->fixed = 1;
		environment_set_value_double(v,p->getValue());
		return 1;
	}

	return 0;
}

/*****************************************************
 Add a new species to the 
******************************************************/
static struct value *simulation_context_add_species(struct simulation_context *sc, Species *s)
{
	char *name;
	struct value *v;

	if (!(name = strdup(s->getId().c_str())))
	{
		fprintf(stderr,"Not enough memory!\n");
		return NULL;
	}

	if (!(v = environment_add_value(&sc->global_env, name)))
		return NULL;
	
	if (s->isSetInitialAmount())
	{
		v->value = s->getInitialAmount();
		v->molecules = s->getInitialAmount();
	} else
	{
		v->value = s->getInitialConcentration();
	}
	v->fixed = s->getBoundaryCondition();
	return v;
}

/***********************************************/

/*****************************************************
 Gets an AST of the stoichiometry factor given
 species reference.
******************************************************/
static ASTNode *get_stoichiometry_ast(const SpeciesReference *ref)
{
	const struct StoichiometryMath *stoichMath = ref->getStoichiometryMath();
	ASTNode *stoich;

	if (stoichMath != NULL)
	{
		return stoichMath->getMath()->deepCopy();
	}

	if (!(stoich = new ASTNode(AST_REAL)))
		return NULL;

	stoich->setValue(ref->getStoichiometry());

	return stoich;
}

/*****************************************************
 For the given SpeciesReference add the formula
 to the right part of its ODE.
******************************************************/
static void simulation_context_add_reference(struct simulation_context *sc, SpeciesReference *ref, const ASTNode *formula, ASTNodeType_t type)
{
	struct value *v;

	if ((v = environment_get_value(&sc->global_env, ref->getSpecies().c_str())))
	{
		ASTNode *prev, *minus, *copy, *stoich;
		
		if (!(prev = v->node))
		{
			if (!(prev = new ASTNode(AST_REAL)))
				goto nomem;

			prev->setValue(0.0);
		}

		if (!(minus = new ASTNode(type)))
			goto nomem;

		/* Extract stoichiometry factor */
		if (!(stoich = get_stoichiometry_ast(ref)))
			goto nomem;

		if (!(copy = formula->deepCopy()))
			goto nomem;

		/* Build component with stoichiometry fractor (if present) */
		ASTNode *component;
		if (stoich != NULL)
		{
			if (!(component = new ASTNode(AST_TIMES)))
				goto nomem;

			component->addChild(stoich);
			component->addChild(copy);
		} else
		{
			component = copy;
		}

		minus->addChild(prev);
		minus->addChild(component);
		v->node = minus;
	}
	
	return;

	nomem:
		fprintf(stderr,"Not enough memory!");
		exit(-1);
}

/*************************************************
 Replaces the POWER with 
*************************************************/
static void fix_power_function(ASTNode *node)
{
	unsigned int i;

	if (node->getType() == AST_POWER)
		node->setType(AST_FUNCTION_POWER);

	for (i=0;i<node->getNumChildren();i++)
	{
		ASTNode *n = node->getChild(i);
		fix_power_function(n);
	}
}

/*************************************************
 Create a new simulation from an SBML file.

 On failure returns NULL.
*************************************************/
struct simulation_context *simulation_context_create_from_sbml_file(const char *filename)
{
	SBMLParser *parser = NULL;
	SBMLDocument *doc;
	Model *model;

	unsigned int numSpecies;
	unsigned int numReactions;
	unsigned int numParameters;
	unsigned int numEvents;
	unsigned int i,j,k;
	
	struct simulation_context *sc;

	struct value *value_first;

	/* Allocate memory for context */
	if (!(sc = (struct simulation_context*)malloc(sizeof(*sc))))
	{
		fprintf(stderr,"Could not allocate memory\n");
		goto bailout;
	}
	memset(sc,0,sizeof(*sc));

	/* Initialize global environment */
	environment_init(&sc->global_env, NULL);

	/* Prepare SBML */
	if (!(parser = new SBMLParser(filename)))
	{
		fprintf(stderr,"Could not parse \"%s\"\n",filename);
		goto bailout;
	}

	if (!(doc = parser->getSBMLDocument()))
		goto bailout;

	if (verbose)
		parser->debugOutputSBML();

	if (!(model = doc->getModel()))
		return 0;

	/* Important parameters of SBML */
	numSpecies = model->getNumSpecies();
	numReactions = model->getNumReactions();
	numParameters = model->getNumParameters();
	numEvents = model->getNumEvents();

	value_first = NULL;

	/* Gather global parameters */
	for (i=0;i<numParameters;i++)
	{
		Parameter *p = model->getParameter(i);
		environment_add_parameter(&sc->global_env, p);
	}

	/* Gather species */
	for (i=0;i<numSpecies;i++)
	{
		Species *sp = model->getSpecies(i);
		simulation_context_add_species(sc, sp);
	}

	/* Allocate space for reactions */
	sc->num_reactions = numReactions;
	if (!(sc->reactions = (struct reaction*)malloc(sizeof(struct reaction)*sc->num_reactions)))
	{
		fprintf(stderr,"Not enough memory\n");
		goto bailout;
	}
	memset(sc->reactions,0,sizeof(struct reaction)*sc->num_reactions);

	/* Gather parameters of the reactions */
	for (i=0;i<numReactions;i++)
	{
		unsigned int numParameter;
		unsigned int numReactants, numProducts;
		struct reaction *r;

		r = &sc->reactions[i];

		/* Initialize environment of reaction parameters */
//		environment_init(&r->env,&sc->global_env);

		Reaction *reaction = model->getReaction(i);
		const char *reactionName = reaction->getName().c_str();
		KineticLaw *kineticLaw = reaction->getKineticLaw();
		ASTNode *formula = kineticLaw->getMath()->deepCopy();

		fix_power_function(formula);

		numParameter = kineticLaw->getNumParameters();
		numReactants = reaction->getNumReactants();
		numProducts = reaction->getNumProducts();

		sc->reactions[i].num_reactants = numReactants;
		if (!(sc->reactions[i].reactants = (struct reference*)malloc(sizeof(sc->reactions[i].reactants[0])*numReactants)))
			goto bailout;
		sc->reactions[i].num_products = numProducts;
		if (!(sc->reactions[i].products = (struct reference*)malloc(sizeof(sc->reactions[i].reactants[0])*numProducts)))
			goto bailout;

		/* Rename local parameters */
		for (j=0;j<numParameter;j++)
		{
			Parameter *p = kineticLaw->getParameter(j);
			ASTNode *new_param;
			const char *org_name = p->getName().c_str();
			struct value *v;
			char *new_param_name;

			if (!(new_param_name = (char*)malloc(strlen(org_name) + strlen(reactionName) + 20)))
				goto bailout;
			sprintf(new_param_name,"%s_%s_r%d",org_name,reactionName,i);
			if (!(new_param = new ASTNode(AST_NAME)))
				goto bailout;
			new_param->setName(new_param_name);

			if (!(v = environment_add_value(&sc->global_env,new_param_name)))
				goto bailout;
			v->fixed = 1;
			environment_set_value_double(v,p->getValue());

			formula->ReplaceArgument(p->getName(),new_param);
		}

		sc->reactions[i].formula = formula;

		/* Add reactants (decrease the species' concentration */
		for (j=0;j<numReactants;j++)
		{
			struct value *ref_v;

			SpeciesReference *ref = reaction->getReactant(j);
			simulation_context_add_reference(sc, ref, formula, AST_MINUS);
			
			if ((ref_v = environment_get_value(&sc->global_env,ref->getSpecies().c_str())))
			{
				sc->reactions[i].reactants[j].value = ref_v;
				if (!(sc->reactions[i].reactants[j].stoich = get_stoichiometry_ast(ref)))
				{
					fprintf(stderr,"Not enough memory!\n");
					goto bailout;
				}
			}
		}
		
		/* Add products (increase the species' concentration */
		for (j=0;j<numProducts;j++)
		{
			struct value *ref_v;

			SpeciesReference *ref = reaction->getProduct(j);
			simulation_context_add_reference(sc, ref, formula, AST_PLUS);

			if ((ref_v = environment_get_value(&sc->global_env,ref->getSpecies().c_str())))
			{
				sc->reactions[i].products[j].value = ref_v;
				if (!(sc->reactions[i].products[j].stoich = get_stoichiometry_ast(ref)))
				{
					fprintf(stderr,"Not enough memory!\n");
					goto bailout;
				}
			}
		}
	}

	/* Allocate the memory for the events */
	if (!(sc->events = (struct event**)malloc(sizeof(sc->events[0])*numEvents)))
	{
		fprintf(stderr,"Could not allocate memory\n");
		goto bailout;
	}
	sc->num_events = numEvents;
	if (!(sc->events_active = (int*)malloc(sizeof(int)*numEvents)))
	{
		fprintf(stderr,"Could not allocate memory\n");
		goto bailout;
	}

	/* Gather events */
	for (i=0;i<numEvents;i++)
	{
		Event *e = model->getEvent(i);
	 	unsigned int numEventAssignments = e->getNumEventAssignments();
	 	struct event *ev;
	 	
	 	if (!(ev = (struct event*)malloc(sizeof(*ev)+sizeof(struct assignment)*numEventAssignments)))
	 	{
			fprintf(stderr,"Could not allocate memory\n");
			goto bailout;
	 	}

	 	ev->num_assignments = numEventAssignments;
	 	ev->trigger = e->getTrigger()->getMath()->deepCopy();

	 	for (j=0;j<numEventAssignments;j++)
	 	{
	 		EventAssignment *ea = e->getEventAssignment(j);
	 		const char *name = ea->getVariable().c_str();

	 		ev->assignments[j].value = environment_get_value(&sc->global_env, name);
	 			
	 		if (ea->getMath())
	 			ev->assignments[j].math = ea->getMath()->deepCopy();
	 	}
	 	sc->events[i] = ev;
	}

	if (verbose)
	{
		/* Print out ODEs */
		/* TODO: Implement a query mechanism */
		fprintf(stderr,"ODEs:\n");

		for (unsigned int i=0;i<sc->global_env.num_values;i++)
		{
			struct value *v = sc->global_env.values[i];
			if (v->node)
				fprintf(stderr,"%s' = %s\n",v->name,SBML_formulaToString(v->node));
		}
	}

	/* Determine (un)fixed variables and build arrays therefrom */
	for (i=0;i<sc->global_env.num_values;i++)
	{
		struct value *v = sc->global_env.values[i];
		if (!v->node || v->fixed)
			sc->num_fixed++;
	}
	sc->num_unfixed = sc->global_env.num_values - sc->num_fixed;

	if (!(sc->unfixed = (struct value**)malloc(sizeof(sc->unfixed[0])*sc->num_unfixed)))
	{
		fprintf(stderr,"Could not parse \"%s\"\n",filename);
		goto bailout;
	}
	
	if (!(sc->fixed = (struct value **)malloc(sizeof(sc->fixed[0])*sc->num_fixed)))
	{
		fprintf(stderr,"Could not parse \"%s\"\n",filename);
		goto bailout;
	}

	for (i=0,j=0,k=0;i<sc->global_env.num_values;i++)
	{
		struct value *v = sc->global_env.values[i];
		if (!v->node || v->fixed) sc->fixed[k++] = v;
		else sc->unfixed[j++] = v;
	}

	/* Temporary space for dynamic linking stuff */
	if (!(sc->dly = (double*)malloc(sizeof(sc->dly[0])*sc->num_unfixed)))
	{
		fprintf(stderr,"Not enough memory\n");
		exit(-1);
	}
	if (!(sc->dlydot = (double*)malloc(sizeof(sc->dlydot[0])*sc->num_unfixed)))
	{
		fprintf(stderr,"Not enough memory\n");
		exit(-1);
	}

	delete parser;
	return sc;
	
bailout:
	if (parser) delete parser;
	return NULL;
}

static double evaluate(struct environment *sc, const ASTNode *node)
{
//	printf("%p type=%d isOper=%d isNumber=%d\n",node,node->getType(),node->isOperator(),node->isNumber());

	switch (node->getType())
	{
	
		case	AST_PLUS: return evaluate(sc, node->getLeftChild()) + evaluate(sc, node->getRightChild());
		case	AST_MINUS:  return evaluate(sc, node->getLeftChild()) - evaluate(sc, node->getRightChild());
		case	AST_TIMES: return evaluate(sc, node->getLeftChild()) * evaluate(sc, node->getRightChild());
		case	AST_DIVIDE: return evaluate(sc, node->getLeftChild()) / evaluate(sc, node->getRightChild());
		case	AST_FUNCTION_POWER: 
		case	AST_POWER: return pow(evaluate(sc, node->getLeftChild()),evaluate(sc, node->getRightChild()));
		case	AST_INTEGER: return node->getInteger();
		case	AST_REAL: return node->getReal();
//		case	AST_REAL_E: break;
///		case	AST_RATIONAL: break;
		case	AST_NAME: return environment_get_value_by_name(sc, node->getName()); break;
//		case	AST_NAME_TIME: break;
//		case	AST_CONSTANT_E: break;
//		case	AST_CONSTANT_FALSE: break;
//		case	AST_CONSTANT_PI: break;
//		case		AST_CONSTANT_TRUE: 	
//				break;
/*
		case		AST_LAMBDA: 	
				break;
		case		AST_FUNCTION:
				break;
		case		AST_FUNCTION_ABS: 	
				break;
		case		AST_FUNCTION_ARCCOS: 	
				break;
		case		AST_FUNCTION_ARCCOSH: 	
				break;
		case		AST_FUNCTION_ARCCOT:	
				break;
		case		AST_FUNCTION_ARCCOTH: 	
				break;
		case		AST_FUNCTION_ARCCSC:	
				break;
		case		AST_FUNCTION_ARCCSCH: 	
				break;
		case		AST_FUNCTION_ARCSEC:	
				break;
		case		AST_FUNCTION_ARCSECH: 	
				break;
		case		AST_FUNCTION_ARCSIN:	
				break;
		case		AST_FUNCTION_ARCSINH: 	
				break;
		case		AST_FUNCTION_ARCTAN:	
				break;
		case		AST_FUNCTION_ARCTANH:	
				break;
		case		AST_FUNCTION_CEILING: 	
				break;
		case		AST_FUNCTION_COS:	
				break;
		case		AST_FUNCTION_COSH: 	
				break;
		case		AST_FUNCTION_COT:	
				break;
		case		AST_FUNCTION_COTH: 	
				break;
		case		AST_FUNCTION_CSC:	
				break;
		case		AST_FUNCTION_CSCH:
				break;
		case		AST_FUNCTION_DELAY: 	
				break;
		case		AST_FUNCTION_EXP:
				break;
		case		AST_FUNCTION_FACTORIAL: 	
				break;
		case		AST_FUNCTION_FLOOR: 	
				break;
		case AST_FUNCTION_LN:	
				break;
		case		AST_FUNCTION_LOG:
				break;
		case		AST_FUNCTION_PIECEWISE: 	
				break;
		case		AST_FUNCTION_POWER:
				break;
		case		AST_FUNCTION_ROOT: 	
				break;
		case		AST_FUNCTION_SEC:	
				break;
		case		AST_FUNCTION_SECH: 	
				break;
		case		AST_FUNCTION_SIN:	
				break;
		case		AST_FUNCTION_SINH:
				break;
		case		AST_FUNCTION_TAN:	
				break;
		case		AST_FUNCTION_TANH: 	
				break;
		case		AST_LOGICAL_AND:
				break;
		case	AST_LOGICAL_NOT: 	
				break;
		case	AST_LOGICAL_OR: 	
				break;
		case	AST_LOGICAL_XOR: 	
				break;*/
		case	AST_RELATIONAL_EQ:  return evaluate(sc, node->getLeftChild()) == evaluate(sc, node->getRightChild());
		case	AST_RELATIONAL_GEQ: return evaluate(sc, node->getLeftChild()) >= evaluate(sc, node->getRightChild());
		case	AST_RELATIONAL_GT:  return evaluate(sc, node->getLeftChild()) > evaluate(sc, node->getRightChild());
		case	AST_RELATIONAL_LEQ: return evaluate(sc, node->getLeftChild()) <= evaluate(sc, node->getRightChild());
		case	AST_RELATIONAL_LT:  return evaluate(sc, node->getLeftChild()) < evaluate(sc, node->getRightChild());
		case	AST_RELATIONAL_NEQ: return evaluate(sc, node->getLeftChild()) != evaluate(sc, node->getRightChild());

		case	AST_UNKNOWN:
				printf("Unknown\n");
				break;
		default:
				printf("ASTNode of type %d not handled!\n",node->getType());
				break;
	}
	
	return 0;
}

static void print(struct environment *sc, const ASTNode *node)
{
//	printf("%p type=%d isOper=%d isNumber=%d\n",node,node->getType(),node->isOperator(),node->isNumber());

	switch (node->getType())
	{
		case	AST_PLUS: printf("(");print(sc, node->getLeftChild()); printf("+"); print(sc, node->getRightChild()); printf(")");break;
		case	AST_MINUS: print(sc, node->getLeftChild()); printf("-"); print(sc, node->getRightChild()); break;
		case	AST_TIMES: print(sc, node->getLeftChild()); printf("*");print(sc, node->getRightChild()); break;
		case	AST_DIVIDE: print(sc, node->getLeftChild()); printf("/");print(sc, node->getRightChild()); break;
		case	AST_FUNCTION_POWER:
		case	AST_POWER: print(sc, node->getLeftChild()); printf("^"); print(sc, node->getRightChild()); break;
		case	AST_INTEGER: printf("%ld",node->getInteger()); break;
		case	AST_REAL: printf("%g",node->getReal()); break;
		case	AST_NAME:  printf("%g",environment_get_value_by_name(sc, node->getName())); break;

		case	AST_UNKNOWN:
				printf("Unknown\n");
				break;
		default:
				printf("ASTNode of type %d not handled!\n",node->getType());
				break;
	}
}

/*********************************************************
 The right hand side of the ODEs (uses dynamical loaded
 function)
**********************************************************/
static int dlf(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
	unsigned int i;
	struct simulation_context *sc = (struct simulation_context*)f_data;

	/* Update values */
	for (i=0;i<sc->num_unfixed;i++)
	{
		sc->dly[i] = NV_Ith_S(y,i);
	}

	/* Calculate ydot */
	sc->dlrhs(t,sc->dly,sc->dlydot,f_data);

	for (i=0;i<sc->num_unfixed;i++)
		NV_Ith_S(ydot,i) = sc->dlydot[i]; 

	return 0;
}

/*********************************************************
 The right hand side of the ODEs
**********************************************************/
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
	unsigned int i;
	struct value *v;
	struct simulation_context *sc = (struct simulation_context*)f_data;

	/* Update values */
	for (i=0;i<sc->num_unfixed;i++)
	{
		v = sc->unfixed[i];
		v->value = NV_Ith_S(y,i);
	}

	/* Calculate ydot */
	for (i=0;i<sc->num_unfixed;i++)
	{
		v = sc->unfixed[i];
		if (v->node && !v->fixed)
			NV_Ith_S(ydot,i) = evaluate(&sc->global_env, v->node);
		else
			NV_Ith_S(ydot,i) = 0;
/*		
		printf("d%s=%g ",v->name,NV_Ith_S(ydot,i));
		printf("%s=%g ",v->name,v->value);
		if (v->node) print(sc,v->node);
		printf(" (t=%g)\n",t);*/
	}
/*
	for (i=0;i<sc->num_values;i++)
	{
		v = sc->values[i];
		printf("%s=%g ",v->name,v->value);
		if (v->node) print(sc,v->node);
		printf(" (t=%g)\n",t);
	}
	*/
	return 0;
}

/**********************************************************
 Returns an NULL-terminated array with names of the
 variables. Memory is freed upon simulation_context_free()
 call. 
***********************************************************/
char **simulation_get_value_names(struct simulation_context *sc)
{
	unsigned int i;

	if (sc->names) return sc->names;
	
	if (!(sc->names = (char**)malloc(sizeof(char*)*(sc->global_env.num_values+1))))
		return NULL;

	for (i=0;i<sc->global_env.num_values;i++)
		sc->names[i] = sc->global_env.values[i]->name;
	sc->names[i] = NULL;
	return sc->names;
}

/**********************************************************
 Compiles the rhs function.
***********************************************************/
static int simulation_context_prepare_jit(struct simulation_context *sc)
{
	unsigned int i;
	const char *filename = "test.c";
	char *command;
	int rc;
	FILE *out;

	if (!(out = fopen(filename,"w")))
	{
		fprintf(stderr,"Unable to open \"%s\" for output.\n",filename);
		return 0;
	}

	/* Build the source code */
	fprintf(out,"#include <stdio.h>\n");
	fprintf(out,"#include <math.h>\n\n");

	if (sc->num_events == 0)
	{
		fprintf(out,"int rhs(double t, double *y, double *ydot, void *f_data)\n");
		fprintf(out,"{\n");
	
		for (i=0;i<sc->global_env.num_values;i++)
		{
			fprintf(out,"\tdouble %s = %g;\n",sc->global_env.values[i]->name,sc->global_env.values[i]->value);
		}
	
		for (i=0;i<sc->num_unfixed;i++)
		{
			struct value *v = sc->unfixed[i];
			fprintf(out,"\t%s = y[%d];\n",v->name,i);
		}
		fprintf(out,"\n");
		for (i=0;i<sc->num_unfixed;i++)
		{
			struct value *v = sc->unfixed[i];
			char *formula = SBML_formulaToString(v->node);
			fprintf(out,"\tydot[%d]=%s;\n",i,formula);
		}
		fprintf(out,"}\n\n");
		fprintf(out,"int check_trigger(double t, double *y, int *events_active)\n{return 0;\n}\n");
	} else
	{
		fprintf(out,"#define geq(a,b) ((a)>=(b))\n");
		fprintf(out,"#define leq(a,b) ((a)<=(b))\n");
		fprintf(out,"\n");
		for (i=0;i<sc->num_fixed;i++)
		{
			struct value *v = sc->fixed[i];
			fprintf(out,"static double %s = %g;\n",v->name,v->value);
		}

		/**** rhs() ****/
		/* The arrays span the the unfixed value space */
		fprintf(out,"\nint rhs(double t, double *y, double *ydot, void *f_data)\n");
		fprintf(out,"{\n");
		for (i=0;i<sc->num_unfixed;i++)
		{
			struct value *v = sc->unfixed[i];
			fprintf(out,"\tdouble %s = y[%d];\n",v->name,i);
		}
		fprintf(out,"\n");
		for (i=0;i<sc->num_unfixed;i++)
		{
			struct value *v = sc->unfixed[i];
			char *formula = SBML_formulaToString(v->node);
			fprintf(out,"\tydot[%d]=%s;\n",i,formula);
		}
		fprintf(out,"}\n\n");

		/**** check_trigger() ****/
		/* The arrays span the whole value space */
		fprintf(out,"int check_trigger(double t, double *y, int *events_active)\n{\n");
		
		for (i=0;i<sc->global_env.num_values;i++)
		{
			if (sc->global_env.values[i]->node)
				fprintf(out,"\tdouble %s=y[%d];\n",sc->global_env.values[i]->name,i);
		}
		fprintf(out,"\n\tint fired = 0;\n\n");
		
		for (i=0;i<sc->num_events;i++)
		{
			struct event *ev;

			ev = sc->events[i];

			fprintf(out, "\tif (%s)\n\t{\n",SBML_formulaToString(ev->trigger));
			fprintf(out, "\t\tif (!events_active[%d])\n",i);
			fprintf(out, "\t\t{\n");
			fprintf(out, "\t\t\tevents_active[%d]=1;\n",i);
			fprintf(out, "\t\t\tfired=1;\n");

			for (unsigned j=0;j<ev->num_assignments;j++)
			{
				struct assignment *a = &ev->assignments[j];

				fprintf(out,"\t\t\ty[%d]=%s=%s;\n",a->value->index,a->value->name,SBML_formulaToString(a->math));
			}

			fprintf(out, "\t\t}\n");
			fprintf(out, "\t} else events_active[%d]=0;\n",i);
		}
		fprintf(out,"\treturn fired;\n}\n");
	}
	
	/* Build the program */
	fclose(out);
	if (!(command = (char*)malloc(500)))
		return 0;

	snprintf(command,500,"gcc -O3 -fPIC -shared test.c -o test.so");
	fprintf(stderr,"%s\n",command);

	if (!(rc = system(command)))
	{
		void *handle = dlopen("./test.so",RTLD_NOW);
		if (handle)
		{
			char *error;

			sc->dlrhs = (int (*)(double t, double *y, double *ydot, void *f_data))dlsym(handle,"rhs");
			if (!(error = dlerror()))
			{
				sc->dlcheck_trigger = (int (*)(double t, double *y, int *events_active))dlsym(handle,"check_trigger");
				if (!(error = dlerror()))
				{
					sc->dlhandle = handle;
					free(command);
					return 1;
				}
			}
			fprintf(stderr,"%s\n",error);
			sc->dlrhs = NULL;
			sc->dlcheck_trigger = NULL;
			dlclose(handle);
		} else
		{
			fprintf(stderr,"dlopen() failed\n%s\n",dlerror());
		}
	}
	free(command);
	return 0;
}

/**********************************************************
 Frees all resources allocated 
***********************************************************/
static void simulation_context_finish_jit(struct simulation_context *sc)
{
	if (sc->dlhandle)
	{
		dlclose(sc->dlhandle);
		sc->dlhandle = NULL;
	}

	sc->dlrhs = NULL;
	sc->dlcheck_trigger = NULL;
}

/**********************************************************
 Perform the event handling. This encompasses the
 evaluation of its trigger. An event is fired, iff
 it transists from false to true.
 
 This function returns 0, if no event has occured.
***********************************************************/
static int simulation_context_check_trigger(struct simulation_context *sc, double t, double *value_space)
{
	if (sc->dlcheck_trigger)
	{
		return sc->dlcheck_trigger(t,value_space,sc->events_active);
	} else
	{
		unsigned int i,j;
		unsigned int fired = 0;
		
		for (i=0;i<sc->num_events;i++)
		{
			double trigger;
			struct event *ev;
	
			ev = sc->events[i];
			
			trigger = evaluate(&sc->global_env, ev->trigger);
	
			if (trigger > 0.0)
			{
				if (!sc->events_active[i])
				{
					/* Event was active before, fire the event by executing its assignments */
					for (j=0;j<ev->num_assignments;j++)
					{
						struct assignment *a = &ev->assignments[j]; 
						if (a)
						{
							if (verbose)
								fprintf(stderr,"Setting %s from %lf to %lf\n",a->value->name,a->value->value,evaluate(&sc->global_env,a->math));
	
							value_space[a->value->index] = a->value->value = evaluate(&sc->global_env,a->math);
							fired = 1;						
						}
						
					}
				}
				sc->events_active[i] = 1;
			} else sc->events_active[i] = 0;
		}

		return fired;
	}
}

/**********************************************************
 Calculates binomial
***********************************************************/
static uint64_t binomial(int N, int K)
{
	int n,k;

	uint64_t binomial[N+1][K+1];

	if (K==0 || K==N) return 1;
	if (K==1 || K==N-1) return N;

    /* base cases */
    for (int k = 1; k <= K; k++) binomial[0][k] = 0;
    for (int n = 0; n <= N; n++) binomial[n][0] = 1;

    /* bottom-up dynamic programming */
	for (n=1;n<=N;n++)
		for (k=1;k<=K;k++)
			binomial[n][k] = binomial[n-1][k-1] + binomial[n-1][k];

	return binomial[N][K];
}

/**********************************************************
 Integrates the simulation using stochastic simulator
***********************************************************/
void simulation_integrate_stochastic(struct simulation_context *sc, struct integration_settings *settings)
{
	double t;
	double tmax = settings->time;

	t = 0;
	
	while (t < tmax)
	{
		double a_all = 0;
		double a_sum = 0;

		/* Calculate propensities */
		for (unsigned int i=0;i<sc->num_reactions;i++)
		{
			struct reaction *r = &sc->reactions[i];
			double h_c = 1.0;

			for (unsigned j=0;j<r->num_reactants;j++)
			{
				struct reference *ref = &r->reactants[j];

				h_c *= binomial(ref->value->molecules,evaluate(&sc->global_env,ref->stoich));
			}

			r->h = h_c;
			r->c = evaluate(&sc->global_env, r->formula);
			r->a = r->h * r->c;

			a_all += r->a;
		}

		double r1 = random()/(double)RAND_MAX;
		double r2 = random()/(double)RAND_MAX;
		
		double tau = (1.0/a_all) * log(1.0/r1);

		for (unsigned int i=0;i<sc->num_reactions;i++)
		{
			struct reaction *r = &sc->reactions[i];

			a_sum += r->a;
			if (a_sum >= r2*a_all)
			{
				for (unsigned j=0;j<r->num_reactants;j++)
					r->reactants[j].value->molecules -= evaluate(&sc->global_env,r->reactants[j].stoich); 
				
				for (unsigned j=0;j<r->num_products;j++)
					r->products[j].value->molecules += evaluate(&sc->global_env,r->products[j].stoich);
				
				break;
			}
		}
		t += tau;

		printf("%lf",t);
		for (unsigned int i=0;i<sc->global_env.num_values;i++)
			printf("\t%d",sc->global_env.values[i]->molecules);
		printf("\n");
	}
}

/**********************************************************
 Integrates the simulation using stochastic simulator
***********************************************************/
void simulation_integrate_stochastic_quick(struct simulation_context *sc, struct integration_settings *settings)
{
	unsigned int i,j;

	double t;
	double tmax = settings->time;

	/* Index indicates whether reaction has been changed */
	int changed[sc->num_reactions];

	/* Contains indices of changed reactions */ 
	unsigned int changed_list[sc->num_reactions];
	
	/* Number of entries in changed_list */
	unsigned int changed_list_size = 0;

	unsigned int stoich_mat_cols = sc->global_env.num_values;
	unsigned int *stoich_mat;
	unsigned int stoich_mat_entries = 0;

	if (!(stoich_mat = (unsigned int*)malloc(sc->num_reactions * stoich_mat_cols * sizeof(unsigned int))))
		return;
	memset(stoich_mat,0,sc->num_reactions * stoich_mat_cols * sizeof(unsigned int));

	/* Build up boolean stoich_mat indicating which species are affected by which reaction */
	for (i=0;i<sc->num_reactions;i++)
	{
		struct reaction *r = &sc->reactions[i];

		for (j=0;j<r->num_reactants;j++)
		{
			stoich_mat[i*stoich_mat_cols + r->reactants[j].value->index] = 1;
			stoich_mat_entries++;
		}

		for (j=0;j<r->num_products;j++)
		{
			stoich_mat[i*stoich_mat_cols + r->products[j].value->index] = -1;
			stoich_mat_entries++;
		}
	}

	/* Build a sparse array containing the indicies of reactions where the species contribute to */
	unsigned int pos;
	int *species_participianting_in_which_reactions_flat;
	int *species_participianting_in_which_reactions[sc->global_env.num_values];

	species_participianting_in_which_reactions_flat = (int*)malloc((stoich_mat_entries + sc->global_env.num_values)*sizeof(int));
	
	pos = 0; /* current position in the flat array */

	for (j=0;j<sc->global_env.num_values;j++)
	{
		unsigned int k;
		int *this_species_participianting_in_which_reactions;

		species_participianting_in_which_reactions[j] = &species_participianting_in_which_reactions_flat[pos];
		this_species_participianting_in_which_reactions = species_participianting_in_which_reactions[j];
		k = 0;

		for (i=0;i<sc->num_reactions;i++)
		{
			if (stoich_mat[i*stoich_mat_cols + j])
				this_species_participianting_in_which_reactions[k++] = i;
		}
		this_species_participianting_in_which_reactions[k] = -1;
		pos += k + 1;
	}

	/* Initially, we flag all reactions as changed */
	for (i=0;i<sc->num_reactions;i++)
		changed_list[i] = i;
	changed_list_size = sc->num_reactions;

	t = 0;

	while (t < tmax)
	{
		double a_all, a_sum;

//		printf("%d: ", changed_list_size);
//		for (i=0;i<changed_list_size;i++)
//		{
//			printf("%d ",changed_list[i]);
//		}
//		printf("\n");

		/* Calculate propensities for all changed-flagged entries */
		for (i=0;i<changed_list_size;i++)
		{
			struct reaction * r = &sc->reactions[changed_list[i]];
			double h_c = 1.0;

			/* unflag */
			changed[changed_list[i]] = 0;

			for (unsigned j=0;j<r->num_reactants;j++)
			{
				struct reference *ref = &r->reactants[j];

				h_c *= binomial(ref->value->molecules,evaluate(&sc->global_env,ref->stoich));
			}

			r->h = h_c;
			r->c = evaluate(&sc->global_env, r->formula);
			r->a = r->h * r->c;
		}

		changed_list_size = 0;
		a_all = 0;

		/* Calculate a_all */
		for (i=0;i<sc->num_reactions;i++)
		{
			struct reaction *r = &sc->reactions[i];
			a_all += r->a;
		}

		double r1 = random()/(double)RAND_MAX;
		double r2 = random()/(double)RAND_MAX;
		
		double tau = (1.0/a_all) * log(1.0/r1);

		a_sum = 0;

		for (i=0;i<sc->num_reactions;i++)
		{
			struct reaction *r = &sc->reactions[i];

			a_sum += r->a;
			if (a_sum >= r2*a_all)
			{
				for (unsigned int j=0;j<r->num_reactants;j++)
				{
					struct value *sp = r->reactants[j].value; 
					sp->molecules -= evaluate(&sc->global_env,r->reactants[j].stoich);

//					printf("reactants: ");
					int *spr = species_participianting_in_which_reactions[sp->index];
					for (unsigned int k=0;spr[k]!=-1;k++)
					{
//						printf("%d ",spr[k]);

						if (!changed[spr[k]])
						{
							changed[spr[k]] = 1;
							changed_list[changed_list_size++] = spr[k]; 
						}
					}
//					printf("\n");
				}
				
				for (unsigned int j=0;j<r->num_products;j++)
				{
					struct value *sp = r->products[j].value;
					sp->molecules += evaluate(&sc->global_env,r->products[j].stoich);

//					printf("products: ");
					int *spr = species_participianting_in_which_reactions[sp->index];
					for (unsigned int k=0;spr[k]!=-1;k++)
					{
//						printf("%d ", spr[k]);

						if (!changed[spr[k]])
						{
							changed[spr[k]] = 1;
							changed_list[changed_list_size++] = spr[k]; 
						}
					}
//					printf("\n");
				}

				break;
			}
		}
		t += tau;

		printf("%lf",t);
		for (unsigned int i=0;i<sc->global_env.num_values;i++)
			printf("\t%d",sc->global_env.values[i]->molecules);
		printf("\n");
	}
	
	free(stoich_mat);
}

/**********************************************************
 Integrates the simulation.
***********************************************************/
void simulation_integrate(struct simulation_context *sc, struct integration_settings *settings)
{
	simulation_integrate_stochastic(sc, settings);
	return;
	
	unsigned i,j;
	double tmax = settings->time;
	unsigned int steps = settings->steps;

	unsigned int num_values = sc->global_env.num_values;
	struct value **values = sc->global_env.values;

	struct value **unfixed = sc->unfixed;
	unsigned int num_unfixed = sc->num_unfixed;

	/* Arrayed version of values. Mainly used for the sample function */
	double *value_space;

	if (verbose)
	{
		fprintf(stderr,"\nInitial values:\n");
		for (i=0;i<num_values;i++)
		{
			fprintf(stderr," %s = %g\n",values[i]->name, values[i]->value);
		}
	}

	if (!(value_space = (double*)malloc(sizeof(double)*num_values)))
	{
		fprintf(stderr,"Not enough memory!\n");
		exit(-1);
	}

	int flag;
	void *cvode_mem = NULL;
	N_Vector initial = NULL, yout = NULL;
	realtype abstol = settings->absolute_error;
	realtype tret;
	int (*rhs)(realtype, _generic_N_Vector*, _generic_N_Vector*, void*);
	
	if (settings->stiff)
		cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	else
		cvode_mem = CVodeCreate(CV_ADAMS,CV_FUNCTIONAL);
	
	if (!cvode_mem)
	{
		fprintf(stderr,"CVodeCreate failed!\n");
		exit(-1);
	}

	if (!(initial = N_VNew_Serial(num_unfixed)))
	{
		fprintf(stderr,"N_VNew_Serial() failed!\n");
		goto out;
	}

	if (!(yout = N_VNew_Serial(num_unfixed)))
	{
		fprintf(stderr,"N_VNew_Serial() failed!\n");
		goto out;
	}

	/* Initialize initial values in value_space (for sample function)
	 * and indicies */
	for (i=0;i<num_values;i++)
		value_space[i] = values[i]->value;

	/* Initialize initial values for ode solver */
	for (i=0;i<num_unfixed;i++)
		NV_Ith_S(initial,i) = unfixed[i]->value;
	
	if (!settings->force_interpreted && simulation_context_prepare_jit(sc))
	{
		if (verbose) fprintf(stderr,"Using compiled right-hand side function\n");
		rhs = dlf;
	} else
	{
		if (verbose) fprintf(stderr,"Using interpreted right-hand side function\n");
		rhs = f;
	}

	flag = CVodeMalloc(cvode_mem, rhs, 0, initial, CV_SS, settings->relative_error, &abstol);
	if (flag < 0)
	{
		fprintf(stderr,"CVodeMalloc failed\n");
		exit(-1);
	}

	CVodeSetFdata(cvode_mem,sc);

	flag = CVDense(cvode_mem,num_unfixed);
	if (flag < 0)
	{
		fprintf(stderr,"CVDense failed\n");
		exit(-1);
	}

	simulation_context_check_trigger(sc,0,value_space);
	if (settings->sample_func)
	{
		if (!settings->sample_func(0.0,num_values,value_space))
			goto out;
	}
	
	CVodeSetMaxNumSteps(cvode_mem,1000000);
	
	for (i=1;i<=steps;i++)
	{
		while ((flag = CVode(cvode_mem,tmax*i/steps,yout,&tret,CV_NORMAL))==CV_TOO_MUCH_WORK);

		if (flag < 0)
		{
			fprintf(stderr,"CVode failed (%d)\n",flag);
			goto out;
		}

		/* Update the value space according to the ode solver */
		for (j=0;j<num_unfixed;j++)
			value_space[unfixed[j]->index] = NV_Ith_S(yout,j);

		if (simulation_context_check_trigger(sc,tret,value_space))
		{
			if (verbose)
				fprintf(stderr,"Reinit at %lf\n",tret);

			/* Update initial values to the current values */
			for (j=0;j<num_unfixed;j++)
				NV_Ith_S(initial,j) = value_space[unfixed[j]->index];

			/* Reinitialize the problem */
			flag = CVodeReInit(cvode_mem, rhs, tret, initial, CV_SS, settings->relative_error, &abstol);
			if (flag < 0)
			{
				fprintf(stderr,"CVodeReInit failed (%d)\n",flag);
				goto out;
			}
		}

		if (settings->sample_func)
		{
			if (!settings->sample_func(tmax*i/steps,num_values,value_space))
				goto out;
		}
	}
	
out:
	simulation_context_finish_jit(sc);

	/* Cleanup */
	if (yout) N_VDestroy_Serial(yout);
	if (initial) N_VDestroy_Serial(initial);
	if (cvode_mem) CVodeFree(&cvode_mem);
	free(value_space);
}

/**********************************************************
 Frees memory associated with a simulation.
***********************************************************/
void simulation_context_free(struct simulation_context *sc)
{
#if 0
	if (sc->values)
	{
		for (unsigned int i=0;i<sc->num_values;i++)
		{
			delete sc->values[i]->node;
			free(sc->values[i]->name);
			free(sc->values[i]);
		}
		free(sc->values);
	}
#endif
	if (sc->events)
	{
		for (unsigned int i=0;i<sc->num_events;i++)
		{
			delete sc->events[i]->trigger;
			for (unsigned j=0;j<sc->events[i]->num_assignments;j++)
				delete sc->events[i]->assignments[j].math;
			free(sc->events[i]);
		}
		free(sc->events);
	}
	free(sc->events_active);

	simulation_context_finish_jit(sc);

	if (sc->names)
		free(sc->names);

	if (sc->unfixed)
		free(sc->unfixed);

	if (sc->fixed)
		free(sc->fixed);
	
	if (sc->dly)
		free(sc->dly);

	if (sc->dlydot)
		free(sc->dlydot);

	if (sc->dlhandle)
		dlclose(sc->dlhandle);

	free(sc);
}

/**********************************************************
 Initializes the settings with some default values. 
***********************************************************/
void integration_settings_init(struct integration_settings *settings)
{
	memset(settings,sizeof(*settings),0);

	settings->time = 10;
	settings->steps = 20;
	settings->absolute_error = 1e-9;
	settings->relative_error = 1e-5;
}
