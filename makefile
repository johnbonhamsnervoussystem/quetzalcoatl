CPP = g++
CFLAGS = -std=c++11 -g
FLAGS = -llapack -lblas 
EFLAG = -I /home/kt27/eigen/eigen-eigen-5a0156e40feb
OBJ = main.o common.o constants.o evalm.o hfrout.o integr.o \
      hfwfn.o project.o qtzio.o solver.o tei.o util.o wigner.o 
OBJTST = test.o common.o constants.o evalm.o hfrout.o hfwfn.o \
      integr.o project.o qtzio.o solver.o tei.o util.o wigner.o

quetzalcoatl : $(OBJ)
	$(CPP) $(CFLAGS) $(FLAGS) $(EFLAG) -o quetzalcoatl $(OBJ) 

all: $(OBJ)

$(OBJ) : %.o: %.cpp
	$(CPP) $(CFLAGS) $(FLAGS) $(EFLAG) -c $< -o $@

test.x : test.o
	$(CPP) $(CFLAGS) $(FLAGS) $(EFLAG) -o test.x $(OBJTST)

test.o : integr.cpp
	$(CPP) $(CFLAGS) $(FLAGS) $(EFLAG) -c integr.cpp -o test.o

