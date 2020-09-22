# Simple Makefile for OpenMP project

all : main.cpp
	g++ main.cpp -fopenmp -o app
clean:
	rm *.o
	rm app
