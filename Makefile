#
#   Makefile for making the executable of program MOPAC
#
#
#    Valid Commands of this makefile
#
#	make linux	Makes the MOPAC file for linux
#	make unix	Makes the MOPAC file for unix 
#	make clean	Clean up disk to minimum config
#
EXE		= ../mopac.exe

mopac.exe:
		@echo "Usage: make <platform>"

linux:
		f2c -w *.f 
		gcc -c -malign-double *.c
		gcc -o $(EXE) *.o -lf2c -lm

unix:
		f77 -w -O *.f -o $(EXE) 

linux2:
		gfortran  -fno-automatic -w -O  *.f -o $(EXE) -fno-range-check
		
linux2debug:
		gfortran  -fno-automatic -w -O -g *.f -o $(EXE) -fno-range-check


clean:
		@mv fdate.c fdate.temp
	 	rm -f *.o
		rm -f *.c
		@mv fdate.temp fdate.c
