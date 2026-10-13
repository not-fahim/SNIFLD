
# Variables

FC = gfortran
FLAGS = -g -fcheck=all -Wall
EXEC = snif


#OBJS = global.o output.o ddsolve.o sweep.o inner.o solver.o version.o parser.o checker.o main.o
OBJS = global.o output.o   version.o parser.o checker.o ldsolve.o sweep.o inner.o solver.o main.o 

# Targets


all: $(EXEC)

# Link Step
$(EXEC): $(OBJS)
	$(FC) $(FLAGS) $(OBJS) -o $(EXEC)


%.o: %.f90
	$(FC) $(FLAGS) -c $< -o $@

#main.o: global.o output.o ddsolve.o sweep.o inner.o solver.o version.o parser.o checker.o
main.o: global.o output.o  version.o parser.o checker.o ldsolve.o sweep.o inner.o solver.o

clean:
	rm -f *.o *.mod 
