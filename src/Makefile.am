bin_PROGRAMS=biozim
biozim_SOURCES=\
	biozim.cpp \
	environment.cpp \
	gnuplot_i.c \
	simulation.cpp \
	sbml/SBMLParser.cpp \
	ode/ODESettings.cpp
biozim_LDADD=@LIBSBML_LIBS@ @SUNDIALS_LIBS@ @LDFLAGS@  

INCLUDES=@LIBSBML_CFLAGS@
