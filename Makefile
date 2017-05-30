# $Id: makefile 16 2007-10-18 12:36:47Z centler $
F77 = gfortran
flag2=-c -O4 -x f77-cpp-input -I. -fPIC 

brns: cows.o 
	$(F77) $(flag1) cows.o -o put 

cows.o: cows.f
	$(F77) $(flag2) cows.f

#brns: cows.o gammln.o GAMMQ.o gcf.o gser.o gammp.o
#	$(F77) $(flag1) cows.o gammln.o GAMMQ.o gcf.o gser.o gammp.o -o put -llapack -lblas
#
#cows.o: cows.f
#	$(F77) $(flag2) cows.f
#gammln.o: gammln.f
#	$(F77) $(flag2) gammln.f
#GAMMQ.o: GAMMQ.f
#	$(F77) $(flag2) GAMMQ.f
#gammp.o: gammp.f
#	$(F77) $(flag2) gammp.f
#gcf.o: gcf.f
#	$(F77) $(flag2) gcf.f
#gser.o: gser.f
	
