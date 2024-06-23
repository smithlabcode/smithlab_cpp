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

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

using std::runtime_error;
using std::string;

MappedRead::MappedRead(const string &line) {
  std::istringstream is;
  is.rdbuf()->pubsetbuf(const_cast<char *>(line.c_str()), line.size());

  string chrom, name, tmp;
  size_t start = 0ul, end = 0ul;
  char strand = '\0';
  double score;
  if (is >> chrom >> start >> tmp) {
    if (find_if(tmp.begin(), tmp.end(),
                [](char c) { return !std::isdigit(c); }) == tmp.end()) {
      end = std::stol(tmp);
      if (!(is >> name >> score >> strand >> seq))
        throw runtime_error("bad line in MappedRead file: " + line);
    }
    else {
      name = tmp;
      if (!(is >> score >> strand >> seq))
        throw runtime_error("bad line in MappedRead file: " + line);
      end = start + seq.length();
    }
    r = GenomicRegion(chrom, start, end, name, score, strand);
    is >> scr;
  }
  else
    throw runtime_error("bad line in MappedRead file: " + line);
}

string MappedRead::tostring() const {
  std::ostringstream oss;
  oss << r; // no chaining for the << of GenomicRegion
  oss << '\t' << seq;
  if (!scr.empty())
    oss << '\t' << scr;
  return oss.str();
}
