/*
 *    Part of SMITHLAB_CPP software
 *
 *    Copyright (C) 2010 University of Southern California and
 *                       Andrew D. Smith
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
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "MappedRead.hpp"
#include "smithlab_utils.hpp"

#include <cassert>
#include <fstream>
#include <algorithm>
#include <stdexcept>

using std::string;
using std::vector;
using std::runtime_error;

std::istream&
operator>>(std::istream& the_stream, MappedRead &mr) {

  string buffer;
  if (getline(the_stream, buffer)) {

    std::istringstream is;
    is.rdbuf()->pubsetbuf(const_cast<char*>(buffer.c_str()), buffer.size());

    string chrom, name;
    size_t start = 0ul, end = 0ul;
    char strand = '\0';
    double score;
    if (is >> chrom >> start >> end >> name >> score >> strand >> mr.seq) {
      mr.r = GenomicRegion(chrom, start, end, name, score, strand);
      is >> mr.scr;
    }
    else throw runtime_error("bad line in MappedRead file: " + buffer);
  }
  return the_stream;
}

std::ostream&
operator<<(std::ostream& the_stream, const MappedRead &mr) {
  the_stream << mr.r << '\t' << mr.seq;
  if (!mr.scr.empty())
    the_stream << '\t' << mr.scr;
  return the_stream;
}
