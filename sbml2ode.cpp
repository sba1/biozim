#include <stdio.h>
#include <stdlib.h>

/* Headers oflibsbml */
#include <sbml/SBMLTypes.h>

/* Headers of sundials */
#include <nvector/nvector_serial.h>
#include <cvode/cvode.h>
#include <cvode/cvode_dense.h>

/* Own headers */
#include "sbml/SBMLParser.h"
#include "ode/ODESettings.h"

/***********************************************/

/* This is our state space */

struct value *value_first;

struct value
{
	/* The name of the value */
	const char *name;
	
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

	v->name = p->getId().c_str();
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

	v->name = s->getId().c_str();
	
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


/********************************************************/

static const char *model_filename;
static int verbose;

/*************************************************
 Displays the program's usage and exits
*************************************************/
static void usage(char *name)
{
	fprintf(stderr, "usage: %s [OPTIONS] SBML-File ...\n"
			"Loads the given SBML-File and performs a simulation run.\n"
			"Specify '-' to read from stdin.\n"
			"\t-h, --help      show this help and quit.\n",
			name);

	exit(1);
}


/*************************************************
 Parse command line args
*************************************************/
static void parse_args(int argc, char *argv[])
{
	int i;
	int filename_given = 0;

	for (i=1;i<argc;i++)
	{
		if (!strcmp(argv[i],"-h") && !strcmp(argv[i],"--help"))
		{
			usage(argv[0]);
			exit(-1);
		} else if (!strcmp(argv[i],"--verbose"))
		{
			verbose = 1;
		} else
		{
			filename_given = 1;
			if (strcmp(argv[i],"-"))
				model_filename = argv[i];
		}
	}
	
	if (!filename_given)
	{
		fprintf(stderr,"No filename has been specifed!\n");
		usage(argv[0]);
		exit(-1);		
	}
}

/******************************************************/

struct simulation_context
{
	unsigned int num_values;
	struct value **values;
};

/*************************************************
 Builds the value map.
*************************************************/
void simulation_context_build_value_map(struct simulation_context *sc)
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
	unsigned int i,j;
	
	struct simulation_context *sc;

	if (!(sc = (struct simulation_context*)malloc(sizeof(*sc))))
	{
		fprintf(stderr,"Could not allocate memory\n");
		goto bailout;
	}

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

//	delete parser;
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
		case	AST_INTEGER: break;
		case	AST_REAL: return node->getReal(); break;
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

/*********************************************************
 The right hand side of the ODEs
**********************************************************/
int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
	unsigned int i;
	struct value *v;
	struct simulation_context *sc = (struct simulation_context*)f_data;

	/* Update values */
	for (i=0;i<sc->num_values;i++)
	{
		v = sc->values[i];
		v->value = NV_Ith_S(y,i);
	}

	/* Calculate ydot */
	for (i=0;i<sc->num_values;i++)
	{
		v = sc->values[i];
		if (v->node && !v->fixed)
			NV_Ith_S(ydot,i) = evaluate(sc, v->node);
		else
			NV_Ith_S(ydot,i) = 0;
	}
	
	return 0;
}


/**********************************************************
 Integrates the simulation.
***********************************************************/
void simulation_integrate(struct simulation_context *sc)
{
	unsigned i,j;
	double tmax = 5;
	unsigned int steps = 10;

	unsigned int num_values = sc->num_values;
	struct value **values = sc->values;

	if (verbose)
	{
		printf("\nInitial values:\n");
		for (i=0;i<num_values;i++)
		{
			printf(" %s = %g\n",sc->values[i]->name, sc->values[i]->value);
		}
	}

	void *cvode_mem;
	int flag;
	realtype *real;
	N_Vector vec;

	if (!(real = (realtype*)malloc(sizeof(*real)*num_values)))
	{
		fprintf(stderr,"Not enough memory!\n");
		exit(-1);
	}

	for (i=0;i<num_values;i++)
		real[i] = values[i]->value;

	if (!(cvode_mem = CVodeCreate(CV_ADAMS,CV_FUNCTIONAL)))
	{
		fprintf(stderr,"CVodeCreate failed!\n");
		exit(-1);
	}

	vec = N_VMake_Serial(num_values,real);

	realtype abstol = 1e-20;

	flag = CVodeMalloc(cvode_mem, f, 0, vec, CV_SS, 1.0e-14, &abstol);
	if (flag < 0)
	{
		fprintf(stderr,"CVodeMalloc failed\n");
		exit(-1);
	}

	/* Set the passed user data */
	CVodeSetFdata(cvode_mem,sc);

	flag = CVDense(cvode_mem,num_values);
	if (flag < 0)
	{
		fprintf(stderr,"CVDense failed\n");
		exit(-1);
	}

	realtype tret;
	N_Vector yout = N_VNew_Serial(num_values);
	
	cout << endl;
	cout << "Results" << endl;

	cout << "Time";

	for (j=0;j<num_values;j++)
	{
		cout << "\t" << values[j]->name;
	}
	cout << endl;
	
	for (i=1;i<=steps;i++)
	{
		flag = CVode(cvode_mem,tmax*i/steps,yout,&tret,CV_NORMAL);
		if (flag < 0)
		{
			fprintf(stderr,"CVode failed\n");
			exit(-1);
		}

		cout << tret;
		for (j=0;j<num_values;j++)
		{
			cout << "\t" << NV_Ith_S(yout,j);
		}
		cout << endl;
	}
}

/**********************************************************
 Frees memory associated with a simulation.
***********************************************************/
void simulation_context_free(struct simulation_context *sc)
{
	free(sc);
}

/**********************************************************
 Main Entry
***********************************************************/
int main(int argc, char **argv)
{
	struct simulation_context *sc;

	parse_args(argc, argv);

	if (!(sc = simulation_context_create_from_sbml_file(model_filename)))
		goto bailout;

	simulation_integrate(sc);

	simulation_context_free(sc);
	

	// Now create settings object for integration of model

/*	ODESettings *settings = new ODESettings();
	cout << *settings;
	// Update some settings
	double time = 1000.0;
	int printstep = 100;
	settings->setTimePointSeries(time, printstep);
	settings->setAbsoluteError(1e-9);
	settings->setRelativeError(1e-4);
	settings->setMaximumSteps(1000);
	cout << *settings;
*/

	return EXIT_SUCCESS;
	
bailout:
	return EXIT_FAILURE;
}


