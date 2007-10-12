#include <stdio.h>
#include <stdlib.h>

/* headers from libsbml */
#include <sbml/SBMLTypes.h>

/* Own headers */
#include "sbml/SBMLParser.h"
#include "ode/ODESettings.h"

/***********************************************/

struct value *value_first;

struct value
{
	const char *name;
	double value;
	const SBase *ref;

	ASTNode *node;

	struct value *next;
};

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
	v->ref = p;
	v->next = value_first;
	value_first = v;
}

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
	v->value = s->getInitialAmount();
	v->ref = s;
	v->next = value_first;
	value_first = v;
}

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

static void value_add_reference(SpeciesReference *ref, const ASTNode *formula, ASTNodeType_t type)
{
	struct value *v;
	if ((v = value_get(ref->getSpecies().c_str())))
	{
		ASTNode *prev, *minus, *copy;
		
		if (!(prev = v->node))
		{
			if (!(prev = new ASTNode(AST_REAL)))
			{
				fprintf(stderr,"Not enough memory!\n");
				exit(-1);
			}
			prev->setValue(0.0);
		}

		if (!(minus = new ASTNode(type)))
		{
			fprintf(stderr,"Not enough memory!\n");
			exit(-1);
		}

		if (!(copy = formula->deepCopy()))
		{
			fprintf(stderr,"Not enough memory!\n");
			exit(-1);
		}
		
		minus->addChild(prev);
		minus->addChild(copy);
		v->node = minus;
	}
}

/***********************************************/


static double evaluate(const ASTNode *node)
{
	printf("%p type=%d isOper=%d isNumber=%d\n",node,node->getType(),node->isOperator(),node->isNumber());

	switch (node->getType())
	{
		case	AST_PLUS: return evaluate(node->getLeftChild()) + evaluate(node->getRightChild());
		case	AST_MINUS:  return evaluate(node->getLeftChild()) - evaluate(node->getRightChild());
		case	AST_TIMES: return evaluate(node->getLeftChild()) * evaluate(node->getRightChild());
		case	AST_DIVIDE: return evaluate(node->getLeftChild()) / evaluate(node->getRightChild());
//		case	AST_POWER: break;
//		case	AST_INTEGER: break;
//		case	AST_REAL: break;
//		case	AST_REAL_E: break;
//		case	AST_RATIONAL: break;
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
				printf("default\n");
				break;
	}
	
	return 0;
}

/**********************************************************
 Main Entry
***********************************************************/
int main(void)
{
	const char *modelfile;
	unsigned int numSpecies;
	unsigned int numReactions;
	unsigned int i,j;

	SBMLParser *parser;
	SBMLDocument *doc;
	
	modelfile = "data/basic-model1-forward-l2.xml";
	
	parser = new SBMLParser(modelfile);
	doc = parser->getSBMLDocument();
	parser->debugOutputSBML();

	Model *model;
	
	if (!(model = doc->getModel()))
		return 0;

	numSpecies = model->getNumSpecies();
	numReactions = model->getNumReactions();

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

	for (i=0;i<numSpecies;i++)
	{
		Species *sp = model->getSpecies(i);
		value_add_species(sp);
	}
	
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

#if 0
	struct value *v = value_first;
	while (v)
	{
		if (v->node)
		{
			cout << v->name << ": " << SBML_formulaToString(v->node) << endl;
		}
		v = v->next;
	}
#endif
	
	// Now create settings object for integration of model
	ODESettings *settings = new ODESettings();
	cout << *settings;
	// Update some settings
	double time = 1000.0;
	int printstep = 100;
	settings->setTimePointSeries(time, printstep);
	settings->setAbsoluteError(1e-9);
	settings->setRelativeError(1e-4);
	settings->setMaximumSteps(1000);
	cout << *settings;
	
	return EXIT_SUCCESS;
}


