mpCOMP=/opt/nfs/mpich-3.1/bin/mpicxx
cppCOMP=g++
FLAGS=-Wall -lm
  
EXECUTABLE=RoomsPar.exe
EXECUTABLEM=RoomsParMPE.exe
EXECUTABLES=RoomsSeq.exe
PAR=ParallelRoomPicker.cpp
PARMPE=ParallelRoomPickerMPE.cpp
SEQ=SequentialRoomPicker.cpp

all:
	$(mpCOMP) -DSPRNG_MPI -DUSE_MPI $(FLAGS) $(PAR) -I/opt/nfs/sprng5/include -L/opt/nfs/sprng5/lib -lsprng -o $(EXECUTABLE)

seq:
	$(cppCOMP) $(FLAGS) $(SEQ) -o $(EXECUTABLES)
	
mpe:
	$(mpCOMP) -DSPRNG_MPI -DUSE_MPI  $(FLAGS)  $(PARMPE) -I/opt/nfs/sprng5/include -L/opt/nfs/sprng5/lib -I/opt/nfs/mpe2-2.4.9b/include -L/opt/nfs/mpe2-2.4.9b/lib -lmpe  -lsprng  -o $(EXECUTABLEM)
	

source:
	source ./source_bash.sh

nodes:
	./station_name_list.sh 101 116 > nodes

run:
	/opt/nfs/mpich-3.1/bin/mpiexec -f nodes -np 4 ./$(EXECUTABLE)
	
runmpe:
	/opt/nfs/mpich-3.2/bin/mpiexec -f nodes -np 4 ./$(EXECUTABLEM)


runseq:
	./$(EXECUTABLES)

conv:
	/opt/nfs/mpe2-2.4.9b/bin/clog2TOslog2 mpe-profile.clog2

jump: 
	/opt/nfs/mpe2-2.4.9b/bin/jumpshot mpe-profile.slog2



clean:
	rm *.dat *.exe nodes *.clog2 *.slog2

