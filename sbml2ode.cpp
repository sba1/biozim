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

/* Own headers */
#include "sbml/SBMLParser.h"
#include "ode/ODESettings.h"


int main(void) {
	const char *modelfile = "data/BIOMD0000000010.xml";
	SBMLParser *parser;
	SBMLDocument *doc;
	parser = new SBMLParser(modelfile);
	doc = parser->getSBMLDocument();
	parser->debugOutputSBML();
	
	// Now create settings object for integration of model
	ODESettings *settings = new ODESettings();
	cout << *settings;
	
	return EXIT_SUCCESS;
}


