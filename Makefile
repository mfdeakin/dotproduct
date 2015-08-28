
CXXFLAGS=-g -std=gnu++11 -Wall

dotprod: dotprod.cpp accurate_math.hpp kobbelt.hpp
	${CXX} ${CXXFLAGS} dotprod.cpp -o dotprod
