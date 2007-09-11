//============================================================================
// Name        : sbml2ode.cpp
// Author      : Peter Robinson
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <stdio.h>
#include <stdlib.h>

/* headers from libsbml */
#include <sbml/SBMLTypes.h>

#include "sbml/SBMLParser.h"

int main(void) {
	const char *modelfile = "data/BIOMD0000000010.xml";
	SBMLParser *parser;
	SBMLDocument *doc;
	parser = new SBMLParser(modelfile);
	doc = parser->getSBMLDocument();
	parser->debugOutputSBML();
	
	return EXIT_SUCCESS;
}


