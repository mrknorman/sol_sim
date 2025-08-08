# ---------------------------------------------------------------
# Detect Homebrew’s libomp location (Intel-mac default: /usr/local/opt/libomp)
# ---------------------------------------------------------------
OMP_ROOT := $(shell brew --prefix libomp)

# Compiler, sources, output
CC      := clang
SRC     := ./src/main.c
TARGET  := sim

# ---------------------------------------------------------------
# Optimisation + vectorisation diagnostics
# ---------------------------------------------------------------
OPTFLAGS := -O3 -march=native -ffast-math \
            -Rpass=loop-vectorize \
            -Rpass-missed=loop-vectorize \
            -Rpass=loop-unroll-and-jam \
            -fsave-optimization-record \
            -Wpedantic

# OpenMP flags for Apple-Clang (needs explicit header & runtime)
OMPFLAGS := -Xpreprocessor -fopenmp -fopenmp-simd \
            -I$(OMP_ROOT)/include \
            -L$(OMP_ROOT)/lib -lomp

# Project headers
INCLUDE  := -I./include -I./res

# ---------------------------------------------------------------
# Always rebuild the final executable (overwrites existing file)
# ---------------------------------------------------------------
all: $(TARGET)

$(TARGET): $(SRC) FORCE         # <- FORCE means “always run this recipe”
	$(CC) $(SRC) $(OPTFLAGS) $(OMPFLAGS) $(INCLUDE) -o $@

# Dummy phony target that is never up-to-date
FORCE:

# ---------------------------------------------------------------
# House-keeping
# ---------------------------------------------------------------
clean:
	rm -f $(TARGET) *.opt.yaml *.opt *.ll *.o

.PHONY: all clean FORCE
