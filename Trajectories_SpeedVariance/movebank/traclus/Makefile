# make the printenv command
#
CC = g++
DEBUG = -O3
INCLUDE = -I../include
CFLAGS = -Wall -c $(INCLUDE) $(DEBUG)
LDFLAGS = -lm

OBJS = TraClusterDoc.o MDPoint.o Trajectory.o Cluster.o ClusterGen.o

all: $(OBJS) Main.o Param.o
	$(CC) $(LDFLAGS) -o traclus $(OBJS) Main.o
	\cp traclus ../bin
	$(CC) $(LDFLAGS) -o param $(OBJS) Param.o
	\cp param ../bin

.cpp.o:
	$(CC) $(CFLAGS) -c $<

clean:
	\rm -f traclus param *.o
