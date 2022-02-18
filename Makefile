all: dBg-min.cpp
	g++ dBg-min.cpp -DNDEBUG -o dBg-min
	
debug: dBg-min.cpp
	g++ dBg-min.cpp -g -O0 -o dBg-min
