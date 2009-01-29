PROFILE = # -pg
CC = gcc-4.3
CXX = g++-4.3
CCFLAGS =	-O3 $(PROFILE) -g -Wall -fmessage-length=0 -fstrict-aliasing -Wstrict-aliasing
CXXFLAGS =	$(CCFLAGS)

ODE_OBJS = ode/ODESettings.o
SBML_OBJS = sbml/SBMLParser.o
OBJS = $(SBML_OBJS) $(ODE_OBJS) \
	environment.o \
	gnuplot_i.o \
	sbml2ode.o \
	simulation.o


LIBS = -L /home/sba/local/lib -lsundials_cvode -lsundials_nvecserial -lsbml -lm -lstdc++
INC = -I . -I /home/sba/local/include

TARGET =	sbml2ode

%.o: %.c
	$(CC) $(CCFLAGS) $(INC) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

$(TARGET):	$(OBJS)
	$(CXX) -o $(TARGET) $(PROFILE) $(OBJS) $(LIBS)

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
