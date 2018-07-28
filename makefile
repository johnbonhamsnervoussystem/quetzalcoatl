CPP = g++
CFLAGS = -std=c++11 -g
HFLAGS = -I ~/working/Eigen -I $(LINTP)/include/libint2/
OBJ = main.o basis.o binio.o common.o constants.o evalm.o hfrout.o integr.o \
      hfwfn.o obarasaika.o project.o qtzio.o solver.o tei.o util.o wigner.o 

quetzalcoatl : $(OBJ)
	$(CPP) $(CFLAGS) $(HFLAGS) -o quetzalcoatl $(OBJ) $(LFLAGS) 

all: $(OBJ)

$(OBJ) : %.o: %.cpp
	$(CPP) $(CFLAGS) $(HFLAGS) -c $< -o $@

