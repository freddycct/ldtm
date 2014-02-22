GSL_INCLUDE=/Users/freddy/Unix/gsl/include
GSL_LIB=/Users/freddy/Unix/gsl/lib

ldtm: user.o ldtm.cpp
	g++ -I$(GSL_INCLUDE) -L$(GSL_LIB) -lgsl -lgslcblas -o ldtm ldtm.cpp user.o

user.o: user.cpp user.h
	g++ -I $(GSL_INCLUDE) -c user.cpp

clean:
	rm -f user.o ldtm