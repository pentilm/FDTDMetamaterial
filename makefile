CC = g++ 
CFLAGS = -Wall -g -ansi
CLIBS = -lm 

files = PML.o fcns.o Algrthm.o drive.o


PML:${files} 
	$(CC) $(CFLAGS) ${files} $(CLIBS) -o PML

%.o: %.cpp
	$(CC) $(CFLAGS)  $(CLIBS)-c $< 

%.o: %.f
	gfortran -c $< 


rmo:
	rm *.o

clean:
	rm *.o *.exe*
