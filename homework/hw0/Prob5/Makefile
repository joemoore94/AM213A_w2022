CC = gcc
CFLAGS = -g -Wall

ApproxPI: Main.o
	$(CC) $(CFLAGS) Main.o -lm -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f ApproxPI *.o *.d *~
