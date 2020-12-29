# Cross-platform Makefile for CLI apps by Scott Carlos Carrion
#
# This program is copyrighted and may not be duplicated or distributed
# without express permission by Scott Carlos Carrion
#
# Copyright (c) Scott Carlos Carrion, 2020
#
# scott.carrion@tamu.edu
#


# Compiler Selection
#CXX = g++
CXX = clang++
#CXX = /scratch/user/scott.carrion/clang-offloading/install/bin/clang++

# Build directory selection (do not include / character)
BUILD_DIR = build

# Executable output
EXE = $(BUILD_DIR)/skyview_omp
#EXE_2 = $(BUILD_DIR)/template_app_2  # Add more than one executable like this

# Source Paths (EXE)
SOURCES = main.cpp  # Main
# Do "SOURCES += path/to/my/other/files" to add other files for compilation of EXE

# Source Paths (EXE_2)
#SOURCES_2 = src/main_2.cpp
# Do "SOURCES_2 += path/to/my/other/files" to add other files for compilation of EXE_2

# Object name generation (EXE)
OBJS = $(addsuffix .o, $(basename $(addprefix $(BUILD_DIR)/, $(notdir $(SOURCES)))))

# Object name generation (EXE_2)
#OBJS_2 = $(addsuffix .o, $(basename $(addprefix $(BUILD_DIR)/, $(notdir $(SOURCES_2)))))

# Extra fun stuff

# Phony targets
.PHONY: all clean
UNAME_S := $(shell uname -s)

# Compiler flags and libraries to link go here
CXXFLAGS = -g -Wall -std=c++11 -Xopenmp-target -march=sm_37 -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda
#CXXFLAGS = -g -std=c++11 -fopenmp
LIBS =

####################################################################################

# Build-specific parameters

ifeq ($(UNAME_S), Linux)  # If platform is Linux
         PLATFORM_NAME = "Linux"
         # Other stuff goes here...
endif

ifeq ($(UNAME_S), Darwin)  # If platform is Mac OS X
         PLATFORM_NAME = "Mac OS X"
        # Other stuff goes here...
endif
####################################################################################

# General build rules

$(BUILD_DIR)/%.o:%.cpp
	mkdir -p $(@D)
	@echo Compiling .cpp file
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# $(BUILD_DIR)/%.o:my/source/folder/%.cpp  # This is a template for adding other sources
#       @echo Compiling .cpp file in this cool directory
#       $(CXX) $(CXXFLAGS) -c -o $@ $<

all: $(EXE) # $(EXE_2)
	@echo Build complete for $(PLATFORM_NAME)

$(EXE): $(OBJS)
	@echo Linking object files \($(EXE)\)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

#$(EXE_2): $(OBJS_2)
#       @echo Linking object files \($(EXE_2)\)
#       $(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

clean:
	rm -rf $(BUILD_DIR)
	-rm skyview*
	-rm prominence*
