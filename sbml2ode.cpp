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
 Return the actual value of the value given by the
 name.
******************************************************/
static double value_get_value(const char *name)
{
	struct value *p = value_first;
	while (p)
	{
		if (!strcmp(name,p->name))
			return p->value;
		p = p->next;
	}
	fprintf(stderr,"value %s not found\n",name);
	return 0;
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


static double evaluate(const ASTNode *node)
{
//	printf("%p type=%d isOper=%d isNumber=%d\n",node,node->getType(),node->isOperator(),node->isNumber());

	switch (node->getType())
	{
	
		case	AST_PLUS: return evaluate(node->getLeftChild()) + evaluate(node->getRightChild());
		case	AST_MINUS:  return evaluate(node->getLeftChild()) - evaluate(node->getRightChild());
		case	AST_TIMES: return evaluate(node->getLeftChild()) * evaluate(node->getRightChild());
		case	AST_DIVIDE: return evaluate(node->getLeftChild()) / evaluate(node->getRightChild());
		case	AST_FUNCTION_POWER:
		case	AST_POWER: return pow(evaluate(node->getLeftChild()),evaluate(node->getRightChild()));
		case	AST_INTEGER: break;
		case	AST_REAL: return node->getReal(); break;
//		case	AST_REAL_E: break;
///		case	AST_RATIONAL: break;
		case	AST_NAME: return value_get_value(node->getName()); break;
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

	/* Update values */
	i = 0;
	v = value_first;
	while (v)
	{
		v->value = NV_Ith_S(y,i);
		v = v->next;
		i++;
	}

	/* Calculate ydot */
	i = 0;
	v = value_first;
	while (v)
	{
		if (v->node && !v->fixed)
			NV_Ith_S(ydot,i) = evaluate(v->node);
		else
			NV_Ith_S(ydot,i) = 0;

		v = v->next;
		i++;
	}
	
	return 0;
}



/**********************************************************
 Main Entry
***********************************************************/
int main(int argc, char **argv)
{
	const char *modelfile;
	unsigned int numSpecies;
	unsigned int numReactions;
	unsigned int i,j;

	SBMLParser *parser;
	SBMLDocument *doc;
	
	if (argc < 2)
	{
		fprintf(stderr,"No SBML-file specified!\n");
		exit(-1);
	}

	modelfile = argv[1];
	
	parser = new SBMLParser(modelfile);
	doc = parser->getSBMLDocument();
	parser->debugOutputSBML();

	Model *model;
	
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

	/* Print out ODEs */
	cout << "ODEs:" << endl;
	unsigned int num_values = 0;
	struct value *v = value_first;
	while (v)
	{
		if (v->node)
			cout << v->name << "' = " << SBML_formulaToString(v->node) << endl;
		v = v->next;
		num_values++;
	}

	/* Convert to an array for easier access */
	struct value **values = (struct value**)malloc(sizeof(*values)*num_values);
	if (!values)
	{
		fprintf(stderr,"Not enough memory\n");
		exit(-1);
	}
	
	v = value_first;
	i = 0;
	while (v)
	{
		values[i++] = v;
		v = v->next;
	}
	
	
	
	{
		double tmax = 20;
		unsigned int steps = 1000;

		cout << endl;
		cout << "Intial values:" << endl;
		for (i=0;i<num_values;i++)
		{
			cout << values[i]->name << "=" << values[i]->value << endl;
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
}


