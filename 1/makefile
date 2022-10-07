#makefile

GCFLAGS=-O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Wcast-qual -Wno-suggest-attribute=format

all: a.out

a.out: gauss.o matrix.o 
	g++ $^

%.o: %.cpp matrix.h
	g++ -c $(GCFLAGS) $< -o $@

clean:
	rm *.o*

rm_out:
	rm *.o
