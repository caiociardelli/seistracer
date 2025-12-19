# Makefile for PyTracer

# Compiler and linker
CC := gcc
# Flags for compiler to build executable
FLAGS1 := \
	-Wall \
	-Wextra \
	-std=c99 \
	-O3
# Flags for compiler to build shared library
FLAGS2 := \
	-std=c99 \
	-Wall \
	-Wextra \
	-shared \
	-fPIC \
	-O3
# Source directory
SRC := src
# Shared objects directory
SHARED := $(SRC)/shared
# Configuration files directory
SETUP := setup
# Configuration files
CONFIG := $(wildcard $(SETUP)/*.h)
# Includes
INC := -I$(SETUP) -I$(SRC)/headers
# Objects directory
OBJ := obj
# Libraries directory
LIB := lib
# Binaries directory
BIN := bin
# Libraries links
LIB_LINK := -lm

# First default target
DEFAULT1 := tables

# Second default target
DEFAULT2 := libtracer.so libio.so libcoefficients.so

# First list of objects
OBJECTS1 := \
	permutations.o \
	amplitudes.o \
	ray.o \
	source_and_model.o \
	io.o

# Second list of objects
OBJECTS2 := \
	ray.o \
	source_and_model.o \
	io.o \
	heapsort.o

# Third list of objects
OBJECTS3 := \
	permutations.o \
	amplitudes.o \
	ray.o \
	source_and_model.o

# Complete path for binaries, libraries, and objects
DFT1  := $(patsubst %, $(BIN)/%, $(DEFAULT1))
DFT2  := $(patsubst %, $(LIB)/%, $(DEFAULT2))
DFT3  := $(patsubst %, $(LIB)/%, $(DEFAULT3))
OBJS1 := $(patsubst %.o, $(OBJ)/%.o, $(OBJECTS1))
OBJS2 := $(patsubst %.o, $(OBJ)/%.o, $(OBJECTS2))
OBJS3 := $(patsubst %.o, $(OBJ)/%.o, $(OBJECTS3))

# Command used for cleaning
RM := rm -rf

#
# Compilation and linking
#
all: objDirectory binDirectory libDirectory $(DFT1) $(DFT2)
	@ echo 'Finished building binary!'

$(BIN)/tables: $(OBJS1) $(OBJ)/tables.o
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS1) $(INC) $^ -o $@ $(LIB_LINK)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(LIB)/libtracer.so: $(OBJS2) $(OBJ)/tracer.o
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS2) $(INC) $^ -o $@ $(LIB_LINK)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(LIB)/libio.so: $(OBJS2) $(OBJ)/io.o
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS2) $(INC) $^ -o $@ $(LIB_LINK)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(LIB)/libcoefficients.so: $(OBJS3) $(OBJ)/coefficients.o
	@ echo 'Building binary using $(CC) linker: $@'
	$(CC) $(FLAGS2) $(INC) $^ -o $@ $(LIB_LINK)
	@ echo 'Finished building binary: $@'
	@ echo ' '

$(OBJ)/permutations.o: $(SHARED)/permutations.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS2) $(INC) -c $(SHARED)/permutations.c -o $@
	@ echo ' '

$(OBJ)/amplitudes.o: $(SHARED)/amplitudes.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS2) $(INC) -c $(SHARED)/amplitudes.c -o $@
	@ echo ' '

$(OBJ)/ray.o: $(SHARED)/ray.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS2) $(INC) -c $(SHARED)/ray.c -o $@
	@ echo ' '

$(OBJ)/source_and_model.o: $(SHARED)/source_and_model.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS2) $(INC) -c $(SHARED)/source_and_model.c -o $@
	@ echo ' '

$(OBJ)/io.o: $(SHARED)/io.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS2) $(INC) -c $(SHARED)/io.c -o $@
	@ echo ' '

$(OBJ)/heapsort.o: $(SHARED)/heapsort.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS2) $(INC) -c $(SHARED)/heapsort.c -o $@
	@ echo ' '

$(OBJ)/coefficients.o: $(SHARED)/coefficients.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS2) $(INC) -c $(SHARED)/coefficients.c -o $@
	@ echo ' '

$(OBJ)/tables.o: $(SRC)/tables.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS1) $(INC) -c $(SRC)/tables.c -o $@
	@ echo ' '

$(OBJ)/tracer.o: $(SRC)/tracer.c $(SETUP)
	@ echo 'Building target using $(CC) compiler: $@'
	$(CC) $(FLAGS2) $(INC) -c $(SRC)/tracer.c -o $@
	@ echo ' '

objDirectory:
	@ mkdir -p $(OBJ)

binDirectory:
	@ mkdir -p $(BIN)

libDirectory:
	@ mkdir -p $(LIB)

clean:
	$(RM) $(OBJ)/ $(BIN)/ $(LIB)/

.PHONY: all clean
