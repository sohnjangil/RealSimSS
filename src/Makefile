all:
	g++ -Wall -g -std=c++11 -c sim_SV.cpp -o sim_SV.o
	g++ -Wall -g -std=c++11 -c realsimss.cpp -o realsimss.o
	g++ -Wall -g -std=c++11 -c RealSimSS.cpp -o RealSimSS.o
	g++ -Wall -g -std=c++11 RealSimSS.o realsimss.o sim_SV.o -o RealSimSS
	mv RealSimSS ../bin
