CPP = g++
#CFLAGS = -std=c++11 -Wall -g -O0
#HFLAGS = -I ~/working/Eigen -I  ~/working/Eigen/unsupported -I /usr/local/lib64
#LFLAGS = -lgsl -lgslcblas -lm
#LFLAGS = -lgsl -lgslcblas -lm
OBJ = main.o qtzio.o

#quetzalcoatl : $(OBJ) 
#	$(CPP) $(CFLAGS) $(HFLAGS) -o quetzalcoatl $(OBJ) $(LFLAGS)
#

quetzalcoatl : $(OBJ) 
	$(CPP) -o quetzalcoatl $(OBJ)

all: $(OBJ)

main.o : %.o: %.cpp
	$(CPP) $(CFLAGS) $(HFLAGS) -c $< -o $@ $(LFLAGS)

qtzio.o : qtzio.cpp
	$(CPP) -I $(BOOST_ROOT)/boost/property_tree/ -c $< -o $@

#$(OBJ) : %.o: %.cpp
#	$(CPP) $(CFLAGS) $(HFLAGS) -c $< -o $@ $(LFLAGS)

