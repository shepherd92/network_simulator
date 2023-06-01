# project name
PROJECT := simplicial_complex

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
CXXFLAGS := -std=c++17 -shared -MMD -MP -Wall -Werror -Wextra -Wno-unused-result -Wsign-compare -pthread -pedantic-errors -save-temps=obj -fwrapv -fstack-protector-strong -fPIC
LDFLAGS  := -lstdc++ -m64 -lpthread -ldl -lutil -lm -L/usr/lib/python3.10/config-3.10-x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -fPIC

# directories
ROOT        := ..
MODULE_ROOT := .
BLDDIR      := $(ROOT)/build
INCDIRNAME  := inc
SRCDIR      := $(MODULE_ROOT)/src
EXTERNAL_INC_DIRS = /usr/include/python3.10/

# definition of extenstions
SRCEXT := cpp
# INLEXT := inl
OBJEXT := so

TARGETS := __all release debug clean

# remove files generated when there was an error
.DELETE_ON_ERROR :
# list phony targets
.PHONY : $(TARGETS)

# do print when directory changes
MAKEFLAGS += --no-print-directory

# add .d to make recognized suffixes
SUFFIXES += .d

# find source files
SRCS := $(shell find $(SRCDIR) -type f -iname "*.$(SRCEXT)")
# INLS := $(shell find $(INCDIR) -type f -iname "*.$(INLEXT)")

# create object names from source filenames
EXEOBJS  := $(EXESRCS:%.$(SRCEXT)=$(OBJDIR)/%.$(OBJEXT))

# create dependency file names from object filenames
DEPS := $(EXESRCS:%.$(SRCEXT)=$(OBJDIR)/%.d)

# automatically create dependency files for object files
-include $(DEPS)

# construct include directory compiler flags
INCDIRS := $(ROOT)/$(INCDIRNAME) $(MODULE_ROOT)/$(INCDIRNAME) $(EXTERNAL_INC_DIRS)
CXXFLAGS += $(addprefix -I , $(INCDIRS))

__all : release

# set target specific compiler and linker flags
debug   : CXXFLAGS += -g3 -pg -O0 --coverage -DDBG
debug   : LDFLAGS  += -pg -lgcov
release : CXXFLAGS += -O3

release debug : $(BLDDIR)/$(PROJECT).$(OBJEXT)

$(BLDDIR)/$(PROJECT).$(OBJEXT) : $(SRCS) | $(BLDDIR)
	$(CXX) $(CXXFLAGS) $^ -o $@
	@touch $(BLDDIR)/__init__.py

# make directories
$(BLDDIR):
	@mkdir -p $@

clean :
	@rm -rf *~ core $(BLDDIR)
