# project name
PROJECT := hypergraph

# suppress Makefile output
# .SILENT :

# compiler and linker settings
CXX      := g++
LD       := g++

# std=c++17
# MP:  add phony target for each dependency other than the main file
# MT:  
# MD:  list user and system-defined include files
# MMD: list user-defined include files
# PIC: position independent code = if the code is shared between multiple processes, it is not possible to load the code to the same memory address
# ldl
CXXFLAGS := -std=c++17 -shared -MMD -MP -Wall -Werror -Wextra -Wno-unused-result -Wsign-compare -pthread -pedantic-errors -save-temps=obj -fstack-protector-strong -fPIC
CXXFLAGS += $(shell python3.10 -m pybind11 --includes)
LDFLAGS  := -lstdc++ -m64 -lpthread -ldl -lutil -lm -L/usr/lib/python3.10/config-3.10-x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -fPIC

# directories
ROOT        := ..
MODULE_ROOT := .
BLDDIR      := $(ROOT)/build
INCDIRNAME  := inc
SRCDIR      := $(MODULE_ROOT)/src
TESTDIR     := $(MODULE_ROOT)/test

# definition of extenstions
SRCEXT   := cpp
# INLEXT := inl
OBJEXT   := o
SHOBJEXT := $(shell python3.10-config --extension-suffix)

TARGETS := __all release debug test clean

# remove files generated when there was an error
.DELETE_ON_ERROR :
# list phony targets
.PHONY : $(TARGETS)

# do print when directory changes
MAKEFLAGS += --no-print-directory

# add .d to make recognized suffixes
SUFFIXES += .d

# find source files
SRCS     := $(shell find $(SRCDIR)  -type f -name "*.$(SRCEXT)")
TESTSRCS := $(shell find $(TESTDIR) -type f -name "*.$(SRCEXT)")
# INLS := $(shell find $(INCDIR) -type f -iname "*.$(INLEXT)")

# create dependency file names from object filenames
DEPS := $(EXESRCS:%.$(SRCEXT)=$(OBJDIR)/%.d)

# automatically create dependency files for object files
-include $(DEPS)

# construct include directory compiler flags
INCDIRS := $(ROOT)/$(INCDIRNAME) $(MODULE_ROOT)/$(INCDIRNAME)
CXXFLAGS += $(addprefix -I , $(INCDIRS))

__all : release

# set target specific compiler and linker flags
debug test : CXXFLAGS += -g3 -pg -O0 --coverage -DDBG -lpython3.10
debug test : LDFLAGS  += -pg -lgcov -lpython3.10
release    : CXXFLAGS += -O3

release debug : $(BLDDIR)/$(PROJECT)$(SHOBJEXT)
test          : $(BLDDIR)/$(PROJECT)$(SHOBJEXT) $(BLDDIR)/$(PROJECT)_test

$(BLDDIR)/$(PROJECT)$(SHOBJEXT) : $(SRCS) | $(BLDDIR)
	$(CXX) $(CXXFLAGS) $^ -o $@
	@touch $(BLDDIR)/__init__.py

$(BLDDIR)/$(PROJECT)_test : $(BLDDIR)/$(PROJECT)_test.$(OBJEXT) | $(BLDDIR)
	$(CXX) $(LDFLAGS) -lpython3.10 $(shell python3.10 -m pybind11 --includes) -l adrcm.cpython-310-x86_64-linux-gnu.so -L$(BLDDIR) $^ -o $@

$(BLDDIR)/$(PROJECT)_test.$(OBJEXT) : $(TESTSRCS) | $(BLDDIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# make directories
$(BLDDIR):
	@mkdir -p $@

clean :
	@rm -rf *~ core $(BLDDIR)
