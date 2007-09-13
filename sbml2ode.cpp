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

/*
 * For now, this main program tries to reproduce integrate.c from SBML_ODESolver
 */
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


