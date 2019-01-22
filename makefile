

cc=gcc -fPIC -shared -lm -Wall -O3

all:
	$(cc) fill_array_2d.c -o _lib_fill_array_2d.so

clean:
	rm _lib_fill_array_2d.so
