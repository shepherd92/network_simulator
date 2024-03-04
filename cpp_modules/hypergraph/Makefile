# project name
PROJECT := $(notdir $(shell pwd))

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
# ltbb: use parallel algorithms, i.e. parallel for_each
CXXFLAGS := -std=c++17 -MMD -MP -Wall -Werror -Wextra -Wno-unused-result -Wsign-compare -pthread -pedantic-errors -save-temps=obj -fstack-protector-strong -fPIC
CXXFLAGS += $(shell python$(PYTHON_VERSION) -m pybind11 --includes)
LDFLAGS  := -lstdc++ -shared -m64 -lpthread -lutil -fPIC
LDFLAGS  += $(shell python$(PYTHON_VERSION)-config --ldflags)

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
SHOBJEXT := $(shell python$(PYTHON_VERSION)-config --extension-suffix)

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
# LIBDIRS := /lib/x86_64-linux-gnu/ /usr/lib/x86_64-linux-gnu/
# LDFLAGS += $(addprefix -L , $(LIBDIRS))

all : release

# set target specific compiler and linker flags
debug test : CXXFLAGS += -g3 -pg -O0 --coverage -DDEBUG -lpython$(PYTHON_VERSION)
debug test : LDFLAGS  += -pg -lgcov -lpython$(PYTHON_VERSION) -ltbb
release    : CXXFLAGS += -O3
release	   : LDFLAGS  += -ltbb

release debug : $(BLDDIR)/$(PROJECT)$(SHOBJEXT)
test          : $(BLDDIR)/$(PROJECT)$(SHOBJEXT) $(BLDDIR)/$(PROJECT)_test

$(BLDDIR)/$(PROJECT)$(SHOBJEXT) : $(SRCS) | $(BLDDIR)
	$(CXX) $^ -o $@ $(CXXFLAGS) $(LDFLAGS) 
	@touch $(BLDDIR)/__init__.py

$(BLDDIR)/$(PROJECT)_test : $(BLDDIR)/$(PROJECT)_test.$(OBJEXT) | $(BLDDIR)
	@ln -sf $(BLDDIR)/$(PROJECT)$(SHOBJEXT) $(BLDDIR)/lib$(PROJECT)$(SHOBJEXT)
	$(CXX) $(LDFLAGS) -L $(BLDDIR) -l $(PROJECT)${SHOBJEXT%.*} $^ -o $@

$(BLDDIR)/$(PROJECT)_test.$(OBJEXT) : $(TESTSRCS) | $(BLDDIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# make directories
$(BLDDIR):
	@mkdir -p $@

clean :
	@rm -rf *~ core $(BLDDIR)
