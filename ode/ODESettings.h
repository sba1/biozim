/**
 * ODESettings.h
 * This class manages the settings for the system of ODEs
 * The class is based on the file integratorSettings.h from 
 * SBML_odeSolver v 1.6
 */


#ifndef ODESETTINGS_H_
#define ODESETTINGS_H_

#include <iostream>

using namespace std;

class ODESettings {
public:
	enum ODE_METHOD {ADAMS_MOULTON, BDF};
	enum ITER_METHOD {NEWTON, FUNCTIONAL};
	enum SENSITIVITY_METHOD {SIMULTANEOUS,STAGGERED, STAGGERED1}; 
	enum SENSITIVITY_ANALYSIS {NO_SENSITIVITY_ANALYSIS = 0, USE_SENSITIVITY_ANALYSIS=1};
private:
	double time; /**< Time to which model is integrated */
	int printStep; /**< Number of output steps from 0 to 'time' */
	double *timePoints; /**< Array of points for designated time course (optional) */
	double absTolerance; /**< Absolute tolerance for ODE integration (CVODE library) */
	double relTolerance; /**< Relative tolerance for ODE integration (CVODE library) */
	int maxStepN; /**< Maximum step number of CVode integration */
	enum ODE_METHOD methodOfIntegration; /**< Method of integration for ODE solver */ 
	enum ITER_METHOD methodOfIteration;
	int maxOrder; /**< Maximum order of ADAMS or BDF method */
	enum SENSITIVITY_ANALYSIS sensitivity; /**< If set to 1 (default: 0), use CVODES for sensitivity analysis */
	enum SENSITIVITY_METHOD sensitivityMethod;
	int haltOnEvent; /**< If not 0, halt on event */
	int steadyState; /**< If not 0, stop integration upon an event. */
	int useJacobian; /**< 1: Use Jacobian ASTs or 0: CVODES internal approximation */
	int storeResults; /**< If not 0, store time-course history */
	
public:
	ODESettings();
	~ODESettings();
	friend ostream &operator<<(ostream &out, ODESettings set);
	
private:
	void setDefaults();
	
};


#endif /*ODESETTINGS_H_*/
