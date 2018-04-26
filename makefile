CPP = g++
CFLAGS = -std=c++11
FLAGS = -llapack -lblas
EFLAG = -I /home/kt27/eigen/eigen-eigen-5a0156e40feb
OBJ = main.o common.o evalm.o hfrout.o hfwfn.o qtzio.o solver.o tei.o util.o

quetzalcoatl : $(OBJ)
	$(CPP) $(CFLAGS) $(FLAGS) $(EFLAG) -o quetzalcoatl $(OBJ)

all: $(OBJ)

$(OBJ) : %.o: %.cpp
	$(CPP) $(CFLAGS) $(FLAGS) $(EFLAG) -c $< -o $@

