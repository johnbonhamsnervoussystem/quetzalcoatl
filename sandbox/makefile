CPP = g++
CFLAGS = -std=c++11 -Wall -g -O0
LDFLAGS = -ljsoncpp
JSON_PATH = /home/xiuhtecuhtli/working/jsoncpp/build/debug/lib
EIGEN_PATH = ${EIGEN_ROOT}
OBJ = basis.o

.PHONY: clean

target: qtz_dbg

qtz_dbg: main.o $(OBJ)
	$(CPP) -o qtz_dbg main.o $(OBJ)

main.o: main.cpp
	$(CPP) $(CFLAGS) -I$(EIGEN_PATH) -c $< -o $@

$(OBJ): %.o: %.cpp
	$(CPP) $(CFLAGS) -I$(EIGEN_PATH) -c $< -o $@

clean:
	rm *.o qtz_dbg
