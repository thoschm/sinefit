all:
	g++ -O3 -march=native -o sinefit src/main.cpp -Iinclude -Wall
