/**
 * TimeCourse.h
 * Objects of this class represent the values obtained
 * for some variable over a time course.
 */


#ifndef TIMECOURSE_H_
#define TIMECOURSE_H_


class TimeCourse {
private:
	int ntimepoints; /**> Number of timepoints, including initial conditions */
	char *name;      /**> Name of species being modeled here */
	double *values;  /**> ntimepoint values over the time course */
	double **sensitivity; /**> Sensitivity time courses */
	
public:
	TimeCourse(char *Name);
	TimeCourse(char *Name, int ntimepoints);
	~TimeCourse();
};


#endif /*TIMECOURSE_H_*/
