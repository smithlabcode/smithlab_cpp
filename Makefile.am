# This file is part of smithlab_cpp
#
# Copyright (C) 2010-2019: Andrew D. Smith and Meng Zhou
#
# Author: Andrew D. Smith and Meng Zhou
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

ACLOCAL_AMFLAGS = -I m4

lib_LTLIBRARIES = libsmithlab_cpp.la

libsmithlab_cpp_la_SOURCES = GenomicRegion.cpp MappedRead.cpp		\
OptionParser.cpp QualityScore.cpp RNG.cpp bisulfite_utils.cpp		\
chromosome_utils.cpp sim_utils.cpp smithlab_os.cpp smithlab_utils.cpp	\
zlib_wrapper.cpp

if HAS_HTSLIB
libsmithlab_cpp_la_SOURCES += htslib_wrapper.cpp
endif

include_HEADERS = GenomicRegion.hpp MappedRead.hpp OptionParser.hpp	\
QualityScore.hpp RNG.hpp bisulfite_utils.hpp chromosome_utils.hpp	\
sim_utils.hpp smithlab_os.hpp smithlab_utils.hpp zlib_wrapper.hpp	\
htslib_wrapper.hpp