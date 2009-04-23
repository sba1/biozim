/**
 * SBMLParser.h
 * Header file uses libsbml to input an SBML model.
 */


#ifndef SBMLPARSER_H_
#define SBMLPARSER_H_


#include <sbml/SBMLTypes.h>



class SBMLParser {
	
private:
	SBMLDocument *document;
	
public:
	~SBMLParser();
	SBMLDocument *getSBMLDocument();
	void debugOutputSBML();
	
	
private:
	void inputSBMLDocument(const char *filename);
	void inputSBMLDocumentFromString(const char *str);
	unsigned long getFileSize(const char* filename);
	unsigned long getCurrentMillis();
	
};


#endif /*SBMLPARSER_H_*/
