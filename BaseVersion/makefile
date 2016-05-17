#compile the file
FC := mpiifort 
FCFLAGS = -c -xHost -O2 -fp-model strict -ipo-c -no-prec-div -static
FLFLAGS =
LIBS = tecio64.a -lm -lstdc++

#DIRs
#DIRS = BC Common Geom Input Output Space Time DynaMesh
#source files
SRC = $(wildcard BC/*.F90) $(wildcard Common/*.F90) $(wildcard Geom/*.F90) $(wildcard Input/*.F90) $(wildcard Output/*.F90) $(wildcard Space/*.F90) $(wildcard Time/*.F90) 
OBJ = $(patsubst %.F90, %.o, $(SRC))

#prog name
PROGRAM = TRANS

all: $(PROGRAM)

$(PROGRAM): $(OBJ) main.o Module_Block.o Module_Global.o
	$(FC) $(FLFLAGS) -o $@ $^ $(LIBS)

$(OBJ): block.mod global.mod 

main.o: block.mod global.mod engineVar.mod

block.mod: Module_Block.o

engineVar.mod: Module_Block.o

global.mod: Module_Global.o Module_Block.o block.mod

%.o: %.F90
	$(FC) $(FCFLAGS) -o $@ $<

clean:
	rm -rf *.o */*.o *.mod TRANS

