/**
 * ODESettings.cpp
 * This class manages the settings for the system of ODEs
 * The class is based on the file integratorSettings.h from 
 * SBML_odeSolver v 1.6
 */


#include "ODESettings.h"
#include <iostream>


ODESettings::ODESettings()
{
	setDefaults();
}

ODESettings::~ODESettings()
{
	// Nothing to do!
}


/**<
 * Sets the default parameters for integration. See the .h file for
 * a description of the meaning of the parameters.
 */
void ODESettings::setDefaults()
{
	this->time = 1.0; 
	this->printStep = 10; 
	this->absTolerance = 1e-18; 
	this->relTolerance = 1e-10; 
	this->maxStepN = 10000; 
	this->methodOfIntegration = BDF; 
	this->methodOfIteration = NEWTON;
	this->maxOrder = 5;
	this->sensitivity = NO_SENSITIVITY_ANALYSIS; 
	this->sensitivityMethod = SIMULTANEOUS;
	this->haltOnEvent = 0;
	this->steadyState = 0;  /* i.e., keep integrating even after event */
	this->useJacobian = 1; 
	this->storeResults = 1; 
}

ostream &operator<<(ostream &out, ODESettings set)
{
	out << endl
		<< "INTEGRATION SETTINGS" <<endl
		<< "--Settings for CVODE--" << endl
		<< "Absolute error tolerance: " << set.absTolerance << endl
		<< "Relative error tolerance: " << set.relTolerance << endl
		<< "Maximum number of steps to reach output time "<< set.maxStepN << endl
		<< "Nonlinear solver method: " ;
		if (set.methodOfIntegration == ODESettings::BDF) out << "BDF";
		else if (set.methodOfIntegration == ODESettings::ADAMS_MOULTON) out << "ADAMS MOULTON";
		else out << "Error: Unknown Method";
	out << endl << "Maximum order: " << set.maxOrder << endl
		<< "Iteration method: ";
	if (set.methodOfIteration == ODESettings::NEWTON) out << "Newton";
	else if (set.methodOfIteration == ODESettings::FUNCTIONAL) out << "Functional";
	else out << "Error: Unknown Method";
	out	<< endl;
	return out;
}


/* end of file */
