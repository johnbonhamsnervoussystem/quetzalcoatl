CPP = g++
CFLAGS = -std=c++11 -g
LINTP = /usr/local/libint/2.0.3-stable
LFLAGS = -L $(LINTP)/lib64/ -lint2
HFLAGS = -I ~/working/Eigen -I $(LINTP)/include/libint2/
OBJ = main.o binio.o common.o constants.o evalm.o hfrout.o integr.o \
      hfwfn.o project.o qtzio.o solver.o tei.o util.o wigner.o 

quetzalcoatl : $(OBJ)
	$(CPP) $(CFLAGS) $(HFLAGS) -o quetzalcoatl $(OBJ) $(LFLAGS) 

all: $(OBJ)

$(OBJ) : %.o: %.cpp
	$(CPP) $(CFLAGS) $(HFLAGS) -c $< -o $@

