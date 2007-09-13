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
	this->endtime = 1.0; 
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
	this->timePoints = NULL;
	setTimePointSeries(this->endtime,this->printStep);
}



/**
 * Use this function to overrisde the defaults for
 * endtime and printStep, and to recalculate the timeseries array. These
 * will be the points for which CVODE calculates values.
 */
void ODESettings::setTimePointSeries(double t, int printstep)
{
	int i;
	double *series;
	this->endtime = t;
	this->printStep = printstep;
	if (printstep<2) {
		cerr << "Error [" << __FILE__<<","<<__LINE__
		<<"]: Printstep must have at least 2 time points" << endl;
		exit(-1);
	}
	series = new double[printstep];
	for (i=1; i<=printStep;++i) {
		series[i-1] = i * t/printstep;
	}
	if (this->timePoints != NULL)
		delete this->timePoints;
	this->timePoints = series;
	
}

void ODESettings::setAbsoluteError(double ae){this->absTolerance = ae;}
void ODESettings::setRelativeError(double re){this->relTolerance = re; }
void ODESettings::setMaximumSteps(int maxn) {this->maxStepN = maxn; }

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
	out << "Sensitivity: ";
	if (set.sensitivity == ODESettings::USE_SENSITIVITY_ANALYSIS) {
		out << "yes" <<endl;
		switch (set.sensitivityMethod) {
		case ODESettings::SIMULTANEOUS:
			out << "\tMethod: Simultaneous" << endl;
			break;
		case ODESettings::STAGGERED:
			out << "\tMethod: Staggered" << endl;
			break;
		case ODESettings::STAGGERED1:
					out << "\tMethod: Staggered1" << endl;
					break;
		default:
			out << "Error: Unrecognized method"<<endl;
		}
	} else if (set.sensitivity == ODESettings::NO_SENSITIVITY_ANALYSIS)
		out << "no" <<endl;
	else out << "Error: Unknown Sensitivity Setting"<<endl;
	out <<  "--Settings for sbml2ode--" << endl
		<< "Jacobian: ";
	if (set.useJacobian)
		out << "generate Jacobian" << endl;
	else 
		out << "CVODE\'s internal approximation"<<endl;
	out << "Event handling: ";
	if (set.haltOnEvent)
		out << "Stop integration"<<endl;
	else
		out << "Keep integrating"<<endl;
	out << "Store results: ";
	if (set.storeResults) out << "yes"<< endl;
	else out << "no"<<endl;
	return out;
}


/* end of file */
