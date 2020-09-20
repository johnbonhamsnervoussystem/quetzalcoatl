CPP = g++
CFLAGS = -std=c++11 -Wall -g -O0
LDFLAGS = -L/home/xiuhtecuhtli/working/jsoncpp/build/debug/src/lib_json -ljsoncpp
OBJ = main.o qtzio.o

quetzalcoatl : $(OBJ) 
	$(CPP) -o quetzalcoatl $(OBJ) $(LDFLAGS)

all: $(OBJ)

main.o : %.o: %.cpp
	$(CPP) $(CFLAGS) -L/home/xiuhtecuhtli/working/jsoncpp/build/debug/lib -ljsoncpp -I/home/xiuhtecuhtli/working/jsoncpp/include/json -c $< -o $@ $(LFLAGS)

qtzio.o : qtzio.cpp
	$(CPP) $(CFLAGS) $(LDFLAGS) -I/home/xiuhtecuhtli/working/jsoncpp/include/json -I$(BOOST_ROOT) -c $< -o $@

