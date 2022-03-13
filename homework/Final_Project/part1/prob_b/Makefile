CC = gcc
CFLAGS = -g -Wall

OBJECTS = main.o ColMajorMat.o MatFunc.o

run: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -lm -o $@ -llapack

%.o: ../%.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f run *.o *.d *~
