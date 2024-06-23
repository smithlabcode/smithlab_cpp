# This file is part of smithlab_cpp
#
# Copyright (C) 2010-2024: Andrew D. Smith
#
# Author: Andrew D. Smith
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

EXTRA_DIST = README.md

ACLOCAL_AMFLAGS = -I m4

# For thing we don't want users to override
AM_CXXFLAGS = -Wall -Wextra -Wpedantic

# Users can override this; by default it would get -O2 -g
CXXFLAGS = -O3 -DNDEBUG

lib_LIBRARIES = libsmithlab_cpp.a

libsmithlab_cpp_a_SOURCES = \
	GenomicRegion.cpp \
	MappedRead.cpp \
	OptionParser.cpp \
	QualityScore.cpp \
	chromosome_utils.cpp \
	bisulfite_utils.cpp \
	smithlab_os.cpp \
	smithlab_utils.cpp \
	zlib_wrapper.cpp \
	dna_four_bit.cpp \
	cigar_utils.cpp \
	sam_record.cpp

if ENABLE_HTS
libsmithlab_cpp_a_SOURCES += htslib_wrapper.cpp
endif

include_HEADERS = \
	GenomicRegion.hpp \
	MappedRead.hpp \
	OptionParser.hpp \
	QualityScore.hpp \
	chromosome_utils.hpp \
	bisulfite_utils.hpp \
	smithlab_os.hpp \
	smithlab_utils.hpp \
	zlib_wrapper.hpp \
	dna_four_bit.hpp \
	cigar_utils.hpp \
	sam_record.hpp

if ENABLE_HTS
include_HEADERS += htslib_wrapper.hpp
endif
