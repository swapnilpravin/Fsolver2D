objects = settings.o partitioner.o io.o communicator.o IBM.o schemes.o boundaryConditions.o solver.o channel.o main.o 


FC = $(HOME)/apps/bin/mpif90
FFLAGS = -ffree-line-length-none -O3 -w
# for profiling: -pg 
#
# intel: (for profiling: -pg)
# FFLAGS = -no-wrap-margin -pg -check bounds


channel: $(objects)
	$(FC) $(FFLAGS) -o channel $(objects)

settings.o : settings.f90
	$(FC) $(FFLAGS) -c settings.f90

io.o : io.f90
	$(FC) $(FFLAGS) -c io.f90

partitioner.o : partitioner.f90
	$(FC) $(FFLAGS) -c partitioner.f90

communicator.o : communicator.f90
	$(FC) $(FFLAGS) -c communicator.f90

schemes.o : schemes.f90
	$(FC) $(FFLAGS) -c schemes.f90

boundaryConditions.o : boundaryConditions.f90
	$(FC) $(FFLAGS) -c boundaryConditions.f90 

solver.o : solver.f90
	$(FC) $(FFLAGS) -c solver.f90

IBM.o : IBM.f90
	$(FC) $(FFLAGS) -c IBM.f90

channel.o : channel.f90
	$(FC) $(FFLAGS) -c channel.f90

main.o : main.f90
	$(FC) $(FFLAGS) -c main.f90

# platelet_count.exe: platelet_count.f90
# 	ifort -o platelet_count.exe platelet_count.f90

clean:
	rm *.o *.mod *.exe *.out *.dat *.plt *.tec *~ *pbs.e* *pbs_run.o* *pbs_run.e* core.*
