CC=gcc
CFLAGS=-I.
DEPS = js_cmath.h js_tdoa.h
OBJ = Diamond_TDOA_V1.o js_cmath.o js_tdoa.o

%.o: %.c $ (DEPS)
	$(CC) -c -o $@ $& $< $(CFLAGS)

diamondTDOA: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	$(RM) *.o
	$(RM) js_cmath js_tdoa Diamond_TDOA_V1 diamondTDOA
	$(RM) *~
