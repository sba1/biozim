/* simulation.cpp */

#include <sbml/SBMLTypes.h>

/* Headers of sundials */
#include <nvector/nvector_serial.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

/* Own headers */
#include "sbml/SBMLParser.h"
#include "ode/ODESettings.h"
#include "simulation.h"

/***********************************************/

struct simulation_context
{
	unsigned int num_values;
	struct value **values;

	char **names;

	unsigned int num_unfixed;
	int *unfixed; /* Contains indices to values */
};

extern int verbose;


/***********************************************/

/* This is our state space */

struct value *value_first;

struct value
{
	/* The name of the value */
	char *name;
	
	/* The value's actual value */
	double value;
	
	/* Whether fixed */
	int fixed;
	
	/* The node of ast describing the right part of the ODE */
	ASTNode *node;

	/* The next value */
	struct value *next;
};

/*****************************************************
 Add the given parameter as value
******************************************************/
static void value_add_parameter(Parameter *p)
{
	struct value *v = (struct value*)malloc(sizeof(*v));
	if (!v)
	{
		fprintf(stderr,"Not enough memory!\n");
		exit(-1);
	}
	memset(v,0,sizeof(*v));

	if (!(v->name = strdup(p->getId().c_str())))
	{
		fprintf(stderr,"Not enough memory!\n");
		exit(-1);
	}
	v->value = p->getValue();
	v->next = value_first;
	value_first = v;
}

/*****************************************************
 Add the given species as value
******************************************************/
static void value_add_species(Species *s)
{
	struct value *v = (struct value*)malloc(sizeof(*v));
	if (!v)
	{
		fprintf(stderr,"Not enough memory!\n");
		exit(-1);
	}
	memset(v,0,sizeof(*v));

	if (!(v->name = strdup(s->getId().c_str())))
	{
		fprintf(stderr,"Not enough memory!\n");
		exit(-1);
	}
		
	if (s->isSetInitialAmount())
		v->value = s->getInitialAmount();
	else
		v->value = s->getInitialConcentration();
	v->next = value_first;
	v->fixed = s->getBoundaryCondition();
	value_first = v;
}

/*****************************************************
 Finds the complete value object by the given
 name.
******************************************************/
static struct value *value_get(const char *name)
{
	struct value *p = value_first;
	while (p)
	{
		if (!strcmp(name,p->name))
			return p;
		p = p->next;
	}
	fprintf(stderr,"value %s not found\n",name);
	return NULL;
}

/*****************************************************
 For the given SpeciesReference add the formula
 to the right part of its ODE.
******************************************************/
static void value_add_reference(SpeciesReference *ref, const ASTNode *formula, ASTNodeType_t type)
{
	struct value *v;
	if ((v = value_get(ref->getSpecies().c_str())))
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
		const struct StoichiometryMath *stoichMath = ref->getStoichiometryMath();
		if (stoichMath != NULL)
		{
			if (!(stoich = stoichMath->getMath()->deepCopy()))
				goto nomem;
		} else stoich = NULL;

		if (ref->getStoichiometry() != 1.0)
		{
			if (!(stoich = new ASTNode(AST_REAL)))
				goto nomem;
			stoich->setValue(ref->getStoichiometry());
		}

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

/***********************************************/


/*************************************************
 Builds the value map.
*************************************************/
static void simulation_context_build_value_map(struct simulation_context *sc)
{
	int i, num_values;
	struct value **values;
	struct value *v;

	/* Convert to an array for easier access */

	num_values = 0;
	v = value_first;
	while (v)
	{
		v = v->next;
		num_values++;
	}

	if (!(values = (struct value**)malloc(sizeof(*values)*num_values)))
	{
		fprintf(stderr,"Not enough memory\n");
		exit(-1);
	}

	i = 0;

	v = value_first;
	while (v)
	{
		values[i++] = v;
		v = v->next;
	}
	
	sc->values = values;
	sc->num_values = num_values;
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
	unsigned int i,j;
	
	struct simulation_context *sc;

	if (!(sc = (struct simulation_context*)malloc(sizeof(*sc))))
	{
		fprintf(stderr,"Could not allocate memory\n");
		goto bailout;
	}
	memset(sc,0,sizeof(*sc));

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

	numSpecies = model->getNumSpecies();
	numReactions = model->getNumReactions();
	numParameters = model->getNumParameters();

	/* Gather global parameters */
	for (i=0;i<numParameters;i++)
	{
		Parameter *p = model->getParameter(i);
		value_add_parameter(p);
	}

	/* Gather parameters of the reactions */
	for (i=0;i<numReactions;i++)
	{
		unsigned int numParameter;

		Reaction *reaction = model->getReaction(i);
		KineticLaw *kineticLaw = reaction->getKineticLaw();

		numParameter = kineticLaw->getNumParameters();
		
		for (j=0;j<numParameter;j++)
		{
			Parameter *p = kineticLaw->getParameter(j);
			value_add_parameter(p);
		}
	}

	/* Gather species */
	for (i=0;i<numSpecies;i++)
	{
		Species *sp = model->getSpecies(i);
		value_add_species(sp);
	}

	simulation_context_build_value_map(sc);

	/* Construct ODEs */
	for (i=0;i<numReactions;i++)
	{
		Reaction *reaction = model->getReaction(i);
		KineticLaw *kineticLaw = reaction->getKineticLaw();

		const ASTNode *formula = kineticLaw->getMath();
		
		unsigned int reactants = reaction->getNumReactants();
		unsigned int products = reaction->getNumProducts();

		for (j=0;j<reactants;j++)
		{
			SpeciesReference *ref = reaction->getReactant(j);
			value_add_reference(ref, formula, AST_MINUS);
		}
		
		for (j=0;j<products;j++)
		{
			SpeciesReference *ref = reaction->getProduct(j);
			value_add_reference(ref, formula, AST_PLUS);
		}
	}

	if (verbose)
	{
		/* Print out ODEs */
		printf("ODEs:\n");

		for (unsigned int i=0;i<sc->num_values;i++)
		{
			struct value *v = sc->values[i];
			if (v->node)
				printf("%s' = %s\n",v->name,SBML_formulaToString(v->node));
		}
	}

	/* Find unfixed variables and build a table therefrom */
	for (i=0;i<sc->num_values;i++)
	{
		struct value *v = sc->values[i];
		if (!v->node || v->fixed)
			continue;
		sc->num_unfixed++;
	}

	if (!(sc->unfixed = (int*)malloc(sizeof(int)*sc->num_unfixed)))
	{
		fprintf(stderr,"Could not parse \"%s\"\n",filename);
		goto bailout;
	}

	for (i=0,j=0;i<sc->num_values;i++)
	{
		struct value *v = sc->values[i];
		if (!v->node || v->fixed)
			continue;
		sc->unfixed[j++] = i;
	}

	delete parser;
	return sc;
	
bailout:
	if (parser) delete parser;
	return NULL;
}

static double simulation_context_get_value(struct simulation_context *sc, const char *symbol)
{
	unsigned int i;
	
	for (i=0;i<sc->num_values;i++)
	{
		if (!strcmp(symbol,sc->values[i]->name))
			return sc->values[i]->value;
	}
	fprintf(stderr,"***Warning***: Symbol \"%s\" not found!\n",symbol); 
	return 0;
}

static double evaluate(struct simulation_context *sc, const ASTNode *node)
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
		case	AST_NAME: return simulation_context_get_value(sc, node->getName()); break;
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
				break;
		case	AST_RELATIONAL_EQ: 	
				break;
		case	AST_RELATIONAL_GEQ: 	
				break;
		case	AST_RELATIONAL_GT: 	
				break;
		case	AST_RELATIONAL_LEQ: 	
				break;
		case	AST_RELATIONAL_LT:
				break;
		case	AST_RELATIONAL_NEQ: 	
				break;
*/
		case	AST_UNKNOWN:
				printf("Unknown\n");
				break;
		default:
				printf("ASTNode of type %d not handled!\n",node->getType());
				break;
	}
	
	return 0;
}

static void print(struct simulation_context *sc, const ASTNode *node)
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
		case	AST_NAME:  printf("%g",simulation_context_get_value(sc, node->getName())); break;

		case	AST_UNKNOWN:
				printf("Unknown\n");
				break;
		default:
				printf("ASTNode of type %d not handled!\n",node->getType());
				break;
	}
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
		v = sc->values[sc->unfixed[i]];
		v->value = NV_Ith_S(y,i);
	}

	/* Calculate ydot */
	for (i=0;i<sc->num_unfixed;i++)
	{
		v = sc->values[sc->unfixed[i]];
		if (v->node && !v->fixed)
			NV_Ith_S(ydot,i) = evaluate(sc, v->node);
		else
			NV_Ith_S(ydot,i) = 0;
	}

#if 0
	printf("\n\nt=%g\n",t);
	for (i=0;i<sc->num_values;i++)
	{
		v = sc->values[i];
		printf("%s=%g  ",sc->values[i]->name,NV_Ith_S(ydot,i));
		if (v->node) print(sc,v->node);
		printf("\n");
	}
#endif
	
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
	
	if (!(sc->names = (char**)malloc(sizeof(char*)*(sc->num_values+1))))
		return NULL;

	for (i=0;i<sc->num_values;i++)
		sc->names[i] = sc->values[i]->name;
	sc->names[i] = NULL;
	return sc->names;
}

void simulation_context_compile(struct simulation_context *sc)
{
	FILE *out = stdout;
}


/**********************************************************
 Integrates the simulation.
***********************************************************/
void simulation_integrate(struct simulation_context *sc, struct integration_settings *settings)
{
	unsigned i,j;
	double tmax = settings->time;
	unsigned int steps = settings->steps;

	unsigned int num_values = sc->num_values;
	struct value **values = sc->values;

	int *unfixed = sc->unfixed;
	unsigned int num_unfixed = sc->num_unfixed;

	/* Used for the sample function */
	double *value_space;

	if (verbose)
	{
		printf("\nInitial values:\n");
		for (i=0;i<num_values;i++)
		{
			printf(" %s = %g\n",sc->values[i]->name, sc->values[i]->value);
		}
	}

	if (!(value_space = (double*)malloc(sizeof(double)*sc->num_values)))
	{
		fprintf(stderr,"Not enough memory!\n");
		exit(-1);
	}

	int flag;
	void *cvode_mem = NULL;
	N_Vector initial = NULL, yout = NULL;
	realtype abstol = settings->absolute_error;
	realtype tret;

	if (!(cvode_mem = CVodeCreate(CV_ADAMS,CV_FUNCTIONAL)))
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

	/* Initialize initial values in value_space (for sample function) */
	for (i=0;i<num_values;i++)
		value_space[i] = values[i]->value;

	/* Initialize initial values for ode solver */
	for (i=0;i<num_unfixed;i++)
		NV_Ith_S(initial,i) = value_space[unfixed[i]];

	flag = CVodeMalloc(cvode_mem, f, 0, initial, CV_SS, settings->relative_error, &abstol);
	if (flag < 0)
	{
		fprintf(stderr,"CVodeMalloc failed\n");
		exit(-1);
	}

	simulation_context_compile(sc);

	/* Set the passed user data */
	CVodeSetFdata(cvode_mem,sc);

	flag = CVDense(cvode_mem,num_unfixed);
	if (flag < 0)
	{
		fprintf(stderr,"CVDense failed\n");
		exit(-1);
	}

	if (settings->sample_func)
	{
		if (!settings->sample_func(0.0,num_values,value_space))
			goto out;
	}
	
	CVodeSetMaxNumSteps(cvode_mem,500000);
	
	for (i=1;i<=steps;i++)
	{
		flag = CVode(cvode_mem,tmax*i/steps,yout,&tret,CV_NORMAL);
		if (flag < 0)
		{
			fprintf(stderr,"CVode failed\n");
			goto out;
		}

		if (settings->sample_func)
		{
			/* Fill the value space */
			for (j=0;j<num_unfixed;j++)
				value_space[unfixed[j]] = NV_Ith_S(yout,j);

			if (!settings->sample_func(tmax*i/steps,num_values,value_space))
				goto out;
		}
	}
	
out:
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
	
	if (sc->names)
		free(sc->names);

	if (sc->unfixed)
		free(sc->unfixed);

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
