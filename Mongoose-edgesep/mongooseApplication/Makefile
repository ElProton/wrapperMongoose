CFLAGS = -I../Include -I../SuiteSparse_config -O3
CXXFLAGS = $(CFLAGS)
LDFLAGS = -L../build/lib
LDLIBS = -lmongoose -lm

connectorTest: Mongoose_EdgeCut_Connector.o connectorTest.o mini_spasm.o mmio.o modules.o

.PHONY: clean

clean:
	rm -rf connectorTest *.o