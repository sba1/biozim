/**
 * SBMLParser.h
 * Header file uses libsbml to input an SBML model.
 */


#ifndef SBMLPARSER_H_
#define SBMLPARSER_H_


#include <sbml/SBMLTypes.h>



class SBMLParser {
	
private:
	const char *fname;
	SBMLDocument *document;
	
public:
	SBMLParser(const char *filename);
	SBMLDocument *getSBMLDocument();
	void debugOutputSBML();
	
	
private:
	void inputSBMLDocument();
	unsigned long getFileSize(const char* filename);
	unsigned long getCurrentMillis();
	
};


#endif /*SBMLPARSER_H_*/
