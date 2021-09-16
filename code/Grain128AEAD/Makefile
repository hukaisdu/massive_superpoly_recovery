CC = g++
CFLAGS = -std=c++11 -O3 -lm -lpthread
HEAD = log.h node.h grain.h
OBJ = grain.o main.o
INC = -I/$$GUROBI_HOME/include/
LIB= -L/$$GUROBI_HOME/lib -lgurobi_c++ -lgurobi91 -lm
TARGET = main

all: compile

compile: $(OBJ)
	$(CC) -o $(TARGET) $^ $(CFLAGS) $(INC) $(LIB)

%.o: %.cpp $(HEAD)
	$(CC) -c -o $@ $< $(CFLAGS) $(INC) $(LIB)

clean:
	rm -f *.o $(TARGET) 
