CC = gcc
CFLAGS = -g -Wall

OBJECTS = main.o ColMajorMat.o MatFunc.o

LinAl: $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -lm -o $@

%.o: ../%.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f LinAl *.o *.d *~
