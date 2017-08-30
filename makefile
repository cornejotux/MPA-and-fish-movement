
nodes=1
threads=7

MPIRUN = mpirun## wanna put something good here sometime
##GLFLAGS = -framework OpenGL -framework Cocoa -framework GLUT 

GLFLAGS = -framework OpenGL -framework Cocoa -framework GLUT

MPCC = mpicxx

CC = $(MPCC)


FISHOFILES = fish.o tga.o

##current: tester
##current: initialConditionsTest
current: gridtest
##current: drawtest
##current: gltest
run:

optimization = 

.SUFFIXES: .c .cpp .o

.c.o:
	$(CC) -c $<

.cpp.o:
	$(CC) -c $<

fish.o: fish.cpp fish.h
	$(CC) $(optimization) -c fish.cpp

tga.o: tga.cpp tga.h
	$(CC) $(optimization) -c tga.cpp

mpisimulation: $(FISHOFILES) mpisimulation.cpp
	$(CC) $(optimization) $(FISHOFILES) mpisimulation.cpp -o mpisimulation

gltest: $(FISHOFILES) gltest.cpp
	$(CC) $(optimization) $(GLFLAGS) $(FISHOFILES) gltest.cpp -o gltest -lc

timetrial: $(FISHOFILES) timetrial.cpp
	$(CC) $(optimization) $(FISHOFILES)  timetrial.cpp -o timetrial

drawtest: $(FISHOFILES) drawtest.cpp
	$(CC) $(optimization) $(FISHOFILES)  drawtest.cpp -o drawtest

tester: $(FISHOFILES) tester.cpp
	$(CC) $(optimization) $(FISHOFILES) tester.cpp -o tester

gridtest: $(FISHOFILES) gridtest.cpp
	$(CC) $(optimization) $(FISHOFILES)  gridtest.cpp -o gridtest

initialConditionsTest: $(FISHOFILES) initialConditionsTest.cpp
	$(CC) $(optimization) $(FISHOFILES)  initialConditionsTest.cpp -o initialConditionsTest

gridtestSymmetrical: $(FISHOFILES) gridtestSymmetrical.cpp
	$(CC) $(optimization) $(FISHOFILES)  gridtestSymmetrical.cpp -o gridtestSymmetrical

runworld: world
	$(MPIRUN) ./world -nodes $(nodes) -tasks_per_node $(threads) -rmpool 1 -euilib us -euidevice sn_all -msg_api mpi

clean:
	rm -rf *# *~ *.o

