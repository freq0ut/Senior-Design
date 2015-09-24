INC_DIR=/Users/betio32/Desktop/myPrograms/myC/myLibraries/js_cmath
CC=gcc
CFLAGS=-I$(INC_DIR)
DEPS = $(INC_DIR)/js_cmath.h js_tdoa.h
OBJ = $(INC_DIR)/js_cmath.o js_tdoa.o Diamond_TDOA_V1.o

%.o: %.c $ (DEPS)
	$(CC) -c -o $@ $& $< $(CFLAGS)

diamondTDOA: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)

.PHONY: clean
clean:
	$(RM) *.o $(INC_DIR)/*.o
	$(RM) js_cmath js_tdoa Diamond_TDOA_V1 diamondTDOA
	$(RM) *~ $(INC_DIR)/*~
