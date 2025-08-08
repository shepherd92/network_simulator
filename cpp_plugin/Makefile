# project name is the name of the directory
PROJECT := $(notdir $(shell pwd))
PYTHON_VERSION := 3.12

# suppress Makefile output
# .SILENT :

# set possibe targets
.DEFAULT_GOAL := release
TARGETS := __all release debug test callgrind memcheck clean

# do print when directory changes
MAKEFLAGS += --no-print-directory

# add .d to make recognized suffixes
SUFFIXES += .d

# directories
ROOT        := .
BLDDIR      := $(ROOT)/build
INCDIRNAME  := inc
SRCDIRNAME  := src
TESTDIR     := $(ROOT)/test

# list include directories recursively
# PATH_PARTS := $(subst /, ,$(dir $1))
# PATH_PARTS_WITHOUT_LAST_TWO := $(wordlist 2,$(words $(PATH_PARTS)),x $(PATH_PARTS))
# DIR_PATHS := $(foreach item,$(PATH_PARTS_WITHOUT_LAST_TWO),\
 	$(eval path := $(if $(path),$(path)/)$(item))\
 	$(path))
# INCDIRS := $(addsuffix /inc, $(DIR_PATHS))
define inc_dirs
$(addsuffix /$(INCDIRNAME),\
	$(foreach item,\
		. $(wordlist 2,$(words $(subst /, ,$(dir $1))),x $(subst /, ,$(dir $1))),\
 		$(eval path := $(if $(path),$(path)/)$(item))\
 		$(path)\
	)\
)
$(eval path :=)
endef

###################################################################
###################################################################
ifndef TARGETDIR  # set target directory to debug or release
###################################################################

# set target directory
TARGETDIR := $(BLDDIR)/debug
release callgrind: TARGETDIR := $(BLDDIR)/release

# export variables set to be used by child make processes (TARGETDIR must be exported)
.EXPORT_ALL_VARIABLES :

# for each target run make again (now target directory is set)
$(TARGETS) :
	@${MAKE} $@

###################################################################
else
###################################################################

# compiler and linker settings
CXX      := g++
LD       := g++

# number of threads
MAKEFLAGS += -j8

# std=c++17
# MP:  add phony target for each dependency other than the main file
# MT:  
# MD:  list user and system-defined include files
# MMD: list user-defined include file dependencies
#
# PIC: position independent code = if the code is shared between multiple processes, it is not possible to load the code to the same memory address
# ldl
# ltbb: use parallel algorithms, i.e. parallel for_each
CXXFLAGS := -std=c++20 -MMD -MP -Wall -Werror -Wextra -Wno-unused-result -Wsign-compare -pedantic-errors -save-temps=obj -fstack-protector-strong -fPIC -DBOOST_ALLOW_DEPRECATED_HEADERS
CXXFLAGS += $(shell python$(PYTHON_VERSION) -m pybind11 --includes)
LDFLAGS  := -lstdc++ -shared -m64 -lutil -fPIC -ltbb
LDFLAGS  += $(shell python$(PYTHON_VERSION)-config --ldflags)

# set target specific compiler and linker flags
debug   : CXXFLAGS += -O0 -g3 -DDEBUG
release : CXXFLAGS += -O3

# definition of extenstions
SRCEXT   := cpp
# INLEXT := inl
OBJEXT   := o
SHOBJEXT := $(shell python$(PYTHON_VERSION)-config --extension-suffix)

# remove files generated when there was an error
.DELETE_ON_ERROR :
# list phony targets
.PHONY : $(TARGETS)

# find source files
SOSRCS := $(shell find $(ROOT) -type f -name "*.$(SRCEXT)" -not -path "$(ROOT)/test/*")
SOOBJS := $(SOSRCS:%.$(SRCEXT)=$(TARGETDIR)/%.$(OBJEXT))

TESTSRCS := $(shell find $(ROOT) -type f -name "*.$(SRCEXT)" -not -name "python_module.cpp")
TESTOBJS := $(TESTSRCS:%.$(SRCEXT)=$(TARGETDIR)/%.$(OBJEXT))
# INLS := $(shell find $(INCDIR) -type f -iname "*.$(INLEXT)")

# create dependency file names from object filenames
DEPS := $(EXESRCS:%.$(SRCEXT)=$(OBJDIR)/%.d)

# automatically create dependency files for object files
-include $(DEPS)

# construct include directory compiler flags
INCDIRS := $(shell find $(ROOT) -name $(INCDIRNAME) -type d)
CXXFLAGS += $(addprefix -I , $(INCDIRS))
# LIBDIRS := /lib/x86_64-linux-gnu/ /usr/lib/x86_64-linux-gnu/
# LDFLAGS += $(addprefix -L , $(LIBDIRS))

__all : release

release : $(TARGETDIR)/$(PROJECT)$(SHOBJEXT)
debug   : $(TARGETDIR)/$(PROJECT)_test

test : $(TARGETDIR)/$(PROJECT)_test
	$(TARGETDIR)/$(PROJECT)_test

# link the shared object library
$(TARGETDIR)/$(PROJECT)$(SHOBJEXT) : $(SOOBJS)
	@echo $(LDFLAGS)
	$(CXX) $^ -o $@ $(LDFLAGS)
	@touch $(BLDDIR)/__init__.py $(TARGETDIR)/__init__.py

# run callgrind
callgrind: $(TARGETDIR)/$(PROJECT)_test
	@rm -f $(TARGETDIR)/callgrind.out
	valgrind --tool=callgrind --callgrind-out-file=$(TARGETDIR)/callgrind.out $(TARGETDIR)/$(PROJECT)_test
	callgrind_annotate --tree=both --inclusive=yes --auto=yes --show-percs=yes $(TARGETDIR)/callgrind.out > $(TARGETDIR)/callgrind.txt
	gprof2dot --format=callgrind --output=$(TARGETDIR)/callgrind.dot $(TARGETDIR)/callgrind.out
	callgrind_annotate --inclusive=yes --show-percs=yes $(TARGETDIR)/callgrind.out > $(TARGETDIR)/callgrind.txt
	dot -Tpdf $(TARGETDIR)/callgrind.dot -o $(TARGETDIR)/callgrind.pdf

memcheck: $(TARGETDIR)/$(PROJECT)_test
	valgrind \
		-v \
		--tool=memcheck \
		--log-file=$(TARGETDIR)/valgrind.log \
		--leak-check=full \
		--show-leak-kinds=all \
		--track-origins=yes \
		$(TARGETDIR)/$(PROJECT)_test

# link the test executable
$(TARGETDIR)/$(PROJECT)_test : $(TESTOBJS)
	$(CXX) $^ -o $@ -lstdc++ -m64 -ltbb

# compile source files
# to get include dirs only from parents: $(addprefix -I,$(call inc_dirs,$<))
$(TARGETDIR)/%.$(OBJEXT) : %.$(SRCEXT)
	@mkdir -p $(dir $@)  # make directory
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean :
	@rm -rf *~ core $(BLDDIR)

###################################################################
endif  # TARGETDIR
###################################################################
###################################################################