/* Copyright (C) 2019 Meng Zhou and Andrew D. Smith
 *
 * Author: Meng Zhou and Andrew D. Smith
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#include "zlib_wrapper.hpp"

#include <algorithm>

using std::string;

igzfstream&
getline(igzfstream &in, string &line) {
  if (gzgets(in.fileobj, &in.buf[0], in.chunk_size)) {
      auto eol = std::find(begin(in.buf), end(in.buf), '\n');
      line = string(begin(in.buf), eol);
  }
  return in;
}

igzfstream&
operator>>(igzfstream &in, string &line) {
  return getline(in, line);
}

ogzfstream&
operator<<(ogzfstream &out, const string &line) {
  gzputs(out.fileobj, line.c_str());
  return out;
}

ogzfstream&
operator<<(ogzfstream &out, const char c) {
  gzputc(out.fileobj, c);
  return out;
}

bool
has_gz_ext(const string &filename) {
  const string ext(".gz");
  return filename.size() >= ext.size()
    && filename.compare(filename.size() - ext.size(), ext.size(), ext) == 0;
}
