#CC = /opt/local/bin/mpicc
#ACC = mpicc
#CC = mpicc -dPARALLEL
#ACC = mpicc -DPARALLEL
CC = g++
ACC = gcc

#FLAG = -gdwarf-2
FLAG = -O2 

.C.o: 
	$(CC) $(GSLINC) -o $*.o -c $(FLAG)  $*.C

.c.o:
	$(ACC)	$(COPT) $(CDEFN)	-c $(FLAG)	$*.c

#EXTLIB = -lm -lmpi -lmpi_cxx -lstdc++
EXTLIB = -lm -lstdc++

EXECS = convModel GM1HMM.exe LoLdHMM.exe ChlHMM.exe

all: $(EXECS)

clean:
	rm -f *.o $(EXECS)

#NRMIN = dlinmin.o dbrent.o nrutil.o mnbrak.o f1dim.o df1dim.o nr_lin_min.o linmin.o brent.o frprmn.o dfpmin.o lnsrch.o
NRMIN =

convModel: convModel.o util.o
	$(CC) -o convModel convModel.o util.o

GM1HMM.exe: GM1HMM.o GM1Network.o HMMAlgorithm.o HMMDirective.o util.o HMMForwardBackward.o Fuzzy.o $(NRMIN)
	$(CC) -o GM1HMM.exe GM1HMM.o GM1Network.o HMMAlgorithm.o HMMDirective.o util.o HMMForwardBackward.o Fuzzy.o $(NRMIN) $(LDFLAGS) $(GSL) $(EXTLIB)

LoLdHMM.exe: LoLdHMM.o LoLdNetwork.o HMMAlgorithm.o HMMDirective.o util.o HMMForwardBackward.o Fuzzy.o $(NRMIN)
	$(CC) -o LoLdHMM.exe LoLdHMM.o LoLdNetwork.o HMMAlgorithm.o HMMDirective.o util.o HMMForwardBackward.o Fuzzy.o $(NRMIN) $(LDFLAGS) $(GSL) $(EXTLIB)

ChlHMM.exe: ChlHMM.o ChlNetwork.o HMMAlgorithm.o HMMDirective.o util.o HMMForwardBackward.o Fuzzy.o $(NRMIN)
	$(CC) -o ChlHMM.exe ChlHMM.o ChlNetwork.o HMMAlgorithm.o HMMDirective.o util.o HMMForwardBackward.o Fuzzy.o $(NRMIN) $(LDFLAGS) $(GSL) $(EXTLIB)
