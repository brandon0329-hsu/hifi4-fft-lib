CC = xt-clang
AR = xt-ar
CFLAGS = -I./include -I/home/yu-chen/下載/HiFi4_VFPU_Library_v4_1_0/hifi4_library/include
LDFLAGS = -L. -lAll_fft_NatureDSP_Signal-hifi4_ss_spfpu_7_llvm-Xtensa-release -lm

# Output files
LIB = lib/libfft_wrapper.a
TEST_EXEC = test/test_fft

# Source files
LIB_SRCS = src/fft_wrapper.c
TEST_SRCS = src/main.c

# Object files
LIB_OBJS = $(LIB_SRCS:.c=.o)
TEST_OBJS = $(TEST_SRCS:.c=.o)

.PHONY: all clean test

all: $(LIB) $(TEST_EXEC)

# Build static library
$(LIB): $(LIB_OBJS)
	@mkdir -p lib
	$(AR) rcs $@ $^

# Build test executable
$(TEST_EXEC): $(TEST_OBJS) $(LIB)
	@mkdir -p test
	$(CC) -o $@ $^ $(LDFLAGS)

# Compile source files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: $(TEST_EXEC)
	xt-run $(TEST_EXEC)

clean:
	rm -f $(LIB_OBJS) $(TEST_OBJS) $(LIB) $(TEST_EXEC)