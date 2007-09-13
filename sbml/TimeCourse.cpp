/**
 * TimeCourse.cpp
 * Objects of this class represent the values obtained
 * for some variable over a time course.
 */

#include "TimeCourse.h"

#include <cstdlib>

using namespace std;

TimeCourse::TimeCourse(char *Name)
{
	this->name = Name;
	this->ntimepoints = 0;
	this->values = NULL;
	this->sensitivity = NULL;
}

TimeCourse::TimeCourse(char *Name, int ntimepoints)
{
	this->name = Name;
	this->ntimepoints = ntimepoints;
	this->values = NULL;
	this->sensitivity = NULL;
}

TimeCourse::~TimeCourse()
{
	delete values;
	// sensitivity.
}


/* end of file */
