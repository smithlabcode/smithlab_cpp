# Makefile from smithlab_cpp C++ code library
#
# Copyright (C) 2010-2019 Andrew D. Smith
#
# Authors: Andrew D. Smith
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.

# ADS: this Makefile will be squashed by any created using autotools
# configure script, if run in the directory containing this file.

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp, %.o, $(SOURCES))
HEADERS = $(wildcard *.hpp)
REQUIRES_HTSLIB = htslib_wrapper_deprecated.o htslib_wrapper.o

ifndef HAVE_HTSLIB
NO_HTSLIB := $(filter-out $(REQUIRES_HTSLIB), $(OBJECTS))
OBJECTS = $(NO_HTSLIB)
endif

STATIC_LIB = libsmithlab_cpp.a

CXX = g++
CXXFLAGS = -Wall -std=c++11
OPTFLAGS = -O3
DEBUGFLAGS = -g
override CPPFLAGS += $(INCLUDEARGS)

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
else
CXXFLAGS += $(OPTFLAGS)
endif

all: $(OBJECTS) static

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(CPPFLAGS) $(LDLIBS) $(LDFLAGS)

static: $(OBJECTS)
	ar cr $(STATIC_LIB) $^

clean:
	@-rm -f *.o *.a *~
.PHONY: clean
