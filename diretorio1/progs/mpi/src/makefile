# Source files.
SRC = $(wildcard *.c)

# Executable file.
EXEC = is

# Builds kernel.
all: $(SRC)
	mpicc -g $(SRC) $(CFLAGS) -o $(BINDIR)/$(EXEC).intel $(LIBS)

# Cleans compilation files.
clean:
	rm -f $(BINDIR)/$(EXEC).intel
