# Rules for generic Makefile
#
# CC = compiler name
# DEPS = any additional folders containing header apart from system headers. Default is current directory.
# SRCDIR = directory containing source files. Default is current directory.
# CPPFILES = names of space-separated (without extensions) cpp files to be compiled.
# TARGET= output binary file name.
# LIBS= any libraries to link with (e.g. -lpthread). Default is none.

# debug or release configuration. Default is release.
# 'make dbg=1' will use debug configuration.
ifdef dbg
CFLAGS+=-g3
CONFIGURATION:=debug
else
CFLAGS+=-O3
CONFIGURATION:=release
endif

# Target platform. Default is x86.
# 'make m64bit=1' will target x64 binary.
ifdef m64bit
CFGLAGS+=-m64
PLATFORM:=x64
else
PLATFORM:=x86
endif

BINDEPTH:=.
BINDIR:=$(BINDEPTH)/bin/$(CONFIGURATION)/$(PLATFORM)
ODIR:=./obj/$(CONFIGURATION)/$(PLATFORM)

OBJ=$(patsubst %,$(ODIR)/%.obj,$(CPPFILES))

$(BINDIR)/$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

$(ODIR)/%.obj: $(SRCDIR)/%.cpp
	$(CC) -c -o $@ $< $(DEPS) $(CFLAGS)

# Object and binary files must be placed in obj and bin directories.
# Create bin and obj directories if they do not exist before compiling.
$(OBJ): | $(ODIR) $(BINDIR)

$(ODIR):
	mkdir -p $(ODIR)
$(BINDIR):
	mkdir -p $(BINDIR)

all: $(BINDIR)/$(TARGET)

# Default run.
run: $(BINDIR)/$(TARGET)
	./$(BINDIR)/$(TARGET) $(ARG)

.PHONY: clean

clean:
	rm -rf obj $(BINDEPTH)/bin/debug/x64/$(TARGET) $(BINDEPTH)/bin/debug/x86/$(TARGET) $(BINDEPTH)/bin/release/x64/$(TARGET) $(BINDEPTH)/bin/release/x86/$(TARGET)
