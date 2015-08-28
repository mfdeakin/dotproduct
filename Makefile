
CXXFLAGS=-O3 -std=gnu++11 -Wall

dotprod: dotprod.cpp accurate_math.hpp kobbelt.hpp Makefile
	${CXX} ${CXXFLAGS} dotprod.cpp -o dotprod
