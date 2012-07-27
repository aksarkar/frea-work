CC = gcc
CFLAGS = -O3 -std=c99

all: ranksum kendall_tau

ranksum: ranksum.c
	$(CC) $(CFLAGS) $< -o $@

kendall_tau: kendall_tau.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ranksum
	rm -f kendall_tau
