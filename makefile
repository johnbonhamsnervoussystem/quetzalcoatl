CPP = g++
CFLAGS = -std=c++11 -Wall -g
HFLAGS = -I ~/working/Eigen -I  ~/working/Eigen/unsupported -I /usr/local/lib64
LFLAGS = -lgsl -lgslcblas -lm
OBJ = main.o basis.o binio.o common.o constants.o evalm.o hfrout.o integr.o \
      hfwfn.o obarasaika.o postscf.o project.o qtzcntrl.o qtzio.o wfn.o solver.o tei.o \
      time_dbg.o util.o wigner.o \

quetzalcoatl : $(OBJ) 
	$(CPP) $(CFLAGS) $(HFLAGS) -o quetzalcoatl $(OBJ) $(LFLAGS) 

all: $(OBJ)

$(OBJ) : %.o: %.cpp
	$(CPP) $(CFLAGS) $(HFLAGS) -c $< -o $@ $(LFLAGS)

