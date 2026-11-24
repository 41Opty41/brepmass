CC ?= cc
CFLAGS ?= -std=c99 -Wall -Wextra -Werror -O2 -Iinclude
LDFLAGS ?= -lm

SRCS = src/mass.c src/tessellate.c src/io_mesh.c

all: build/libbrepmass.a build/brepmass build/test_brepmass

build:
	mkdir -p build

build/libbrepmass.a: $(SRCS) | build
	$(CC) $(CFLAGS) -c src/mass.c -o build/mass.o
	$(CC) $(CFLAGS) -c src/tessellate.c -o build/tessellate.o
	$(CC) $(CFLAGS) -c src/io_mesh.c -o build/io_mesh.o
	ar rcs $@ build/mass.o build/tessellate.o build/io_mesh.o

build/brepmass: cli/main.c build/libbrepmass.a | build
	$(CC) $(CFLAGS) cli/main.c build/libbrepmass.a $(LDFLAGS) -o $@

build/test_brepmass: tests/test_brepmass.c build/libbrepmass.a | build
	$(CC) $(CFLAGS) tests/test_brepmass.c build/libbrepmass.a $(LDFLAGS) -o $@

test: build/test_brepmass
	./build/test_brepmass

clean:
	rm -rf build
