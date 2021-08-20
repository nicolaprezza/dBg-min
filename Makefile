all:
	g++ dBg-min.cpp -DNDEBUG -o dBg-min
	
debug:
	g++ dBg-min.cpp -o dBg-min
