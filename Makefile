ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
CXXFLAGS  += $(ROOTCFLAGS)
LIBS       = $(ROOTLIBS) 
GLIBS      = $(ROOTGLIBS)
GXX	   = g++ -Wall -O3


RKdemo: LSQFit.cpp
	$(GXX) LSQFit.cpp -o LSQFit $(ROOTCFLAGS) $(LIBS) $(ROOTGLIBS)

clean:
	rm -f *.d *.so *.pcm LSQFit
