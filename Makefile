CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -fstrict-aliasing -Wstrict-aliasing

ODE_OBJS = ode/ODESettings.o

SBML_OBJS = sbml/SBMLParser.o sbml/TimeCourse.o

OBJS =		sbml2ode.o $(SBML_OBJS) $(ODE_OBJS)

LIBS = -L /home/sba/local/lib -lsundials_cvode -lsundials_nvecserial -lsbml -lm -lstdc++
INC = -I /home/sba/local/include

TARGET =	sbml2ode

%.o: %.cpp
	g++ $(CXXFLAGS) $(INC) -c $< -o $@

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
