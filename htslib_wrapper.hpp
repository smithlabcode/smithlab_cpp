/*    Copyright (C) 2018 Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */

#ifndef HTSLIB_WRAPPER
#define HTSLIB_WRAPPER

#include <string>

#include <htslib/sam.h>

class GenomicRegion;
class MappedRead;

struct HTSFile {
  HTSFile(const std::string &fn);
  ~HTSFile();
  operator bool() const {return is_good;}

  samFile *the_file;
  bam_hdr_t *the_header;
  bool is_good;
};

HTSFile &
operator>>(HTSFile &h, GenomicRegion &r);

HTSFile &
operator>>(HTSFile &h, MappedRead &mr);

#endif
