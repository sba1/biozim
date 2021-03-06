/**
 * SBMLParser.cpp
 * Header file uses libsbml to input an SBML model.
 * Code in this class based on example code from libsbml distribution
 */

#include "SBMLParser.h"
#include <iostream>
#include <stddef.h>
#include <sys/stat.h>
#include <sys/time.h>

extern int verbose;

SBMLParser::~SBMLParser()
{
	if (document != NULL)
		delete document;
}

SBMLDocument *SBMLParser::getSBMLDocument()
{
	return document;
}



void SBMLParser::inputSBMLDocument(const char *fname)
{
	SBMLReader *reader;
	unsigned long start, stop;
	start = getCurrentMillis();
	
	reader = new SBMLReader();
	document = reader->readSBML(fname);
	stop = getCurrentMillis();
	
	unsigned int errors = document->getNumErrors();

	if (verbose)
	{
		cout << endl;
		cout << "            filename: "<< fname            << endl;
		cout << "           file size: "<< getFileSize(fname)     << endl;
		cout << "      read time (ms): " << stop - start          << endl;
		cout << " validation error(s): " << errors                << endl;
		cout << endl;
	}

	document->printErrors(cerr);
	delete reader; 
}

void SBMLParser::inputSBMLDocumentFromString(const char *str)
{
	SBMLReader *reader;
	unsigned long start, stop;
	start = getCurrentMillis();
	
	reader = new SBMLReader();
	document = reader->readSBMLFromString(str);
	stop = getCurrentMillis();
	
	unsigned int errors = document->getNumErrors();

	if (verbose)
	{
		cout << endl;
		cout << "            filename: <stdin>" << endl;
		cout << "           file size: "<< strlen(str)    << endl;
		cout << "      read time (ms): " << stop - start          << endl;
		cout << " validation error(s): " << errors                << endl;
		cout << endl;
	}

	document->printErrors(cerr);
	delete reader; 
	
}

void SBMLParser::debugOutputSBML()
{
	if (this->document == NULL)
		return;

	unsigned int level   = document->getLevel  ();
	unsigned int version = document->getVersion();

	cout << endl
     		<< " (Level " << level << ", version " << version << ")" << endl;

	Model* model = document->getModel();

	if (model == 0)
	{
		cout << "No model present." << endl;
		return;
	}
	cout << "               "
	       << (level == 1 ? "name: " : "  id: ")
	       << (model->isSetId() ? model->getId() : "(empty)") << endl;

	  if (model->isSetSBOTerm())
	    cout << "      model sboTerm: " << model->getSBOTerm() << endl;

	  cout << "functionDefinitions: " << model->getNumFunctionDefinitions() << endl;
	  cout << "    unitDefinitions: " << model->getNumUnitDefinitions    () << endl;
	  cout << "   compartmentTypes: " << model->getNumCompartmentTypes   () << endl;
	  cout << "        specieTypes: " << model->getNumSpeciesTypes       () << endl;
	  cout << "       compartments: " << model->getNumCompartments       () << endl;
	  cout << "            species: " << model->getNumSpecies            () << endl;
	  cout << "         parameters: " << model->getNumParameters         () << endl;
	  cout << " initialAssignments: " << model->getNumInitialAssignments () << endl;
	  cout << "              rules: " << model->getNumRules              () << endl;
	  cout << "        constraints: " << model->getNumConstraints        () << endl;
	  cout << "          reactions: " << model->getNumReactions          () << endl;
	  cout << "             events: " << model->getNumEvents             () << endl;
	  cout << endl;


}

unsigned long SBMLParser::getFileSize(const char* filename)
{
	struct stat s;
	unsigned long result = 0;
	if (stat(filename,&s)==0)
		result = s.st_size;
	return result;
}

unsigned long SBMLParser::getCurrentMillis()
{
	unsigned long result = 0;
	struct timeval tv;
	if (gettimeofday(&tv,NULL)==0)
		result = static_cast<unsigned long>(tv.tv_sec*1000 + tv.tv_usec * .001);
	return result;
}

/* end of file */
