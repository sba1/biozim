CXXFLAGS =	-O2 -g -Wall -fmessage-length=0


ODE_OBJS = ode/ODESettings.o

SBML_OBJS = sbml/SBMLParser.o

OBJS =		sbml2ode.o $(SBML_OBJS) $(ODE_OBJS)

LIBS =	-lsbml -lm -lstdc++

TARGET =	sbml2ode

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
