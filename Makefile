CXXFLAGS =	-O2 -g -Wall -fmessage-length=0

OBJS =		sbml2ode.o sbml/SBMLParser.o

LIBS =	-lsbml -lm -lstdc++

TARGET =	sbml2ode

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
