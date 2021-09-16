CC = g++
CFLAGS = -std=c++11 -O3 -lm -lpthread
DEPS = log.h deg.h node.h trivium.h 
OBJ = deg.o trivium.o main.o
INC = -I/$$GUROBI_HOME/include/
LIB= -L/$$GUROBI_HOME/lib -lgurobi_c++ -lgurobi91 -lm
TARGET = main

all: compile

compile: $(OBJ)
	$(CC) -o $(TARGET) $^ $(CFLAGS) $(INC) $(LIB)

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(INC) $(LIB)

clean:
	/bin/rm -f *.o $(TARGET) 
