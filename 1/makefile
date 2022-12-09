#makefile

CPPFLAGS= -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Wcast-qual -Wno-suggest-attribute=format

.PHONY: all debug

all: a.out

%.o: %.cpp
	g++ -c  $(CPPFLAGS) $< -o $@

a.out: gauss.o matrix.o
	g++ -pthread -o a.out gauss.o matrix.o

clean:
	rm -f *.out *.o
