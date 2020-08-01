/* Copyright (C) 2019 Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * This software is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 */

#include "cigar_utils.hpp"

#include <exception>
#include <sstream>
#include <string>

using std::string;
using std::to_string;
using std::runtime_error;

void
apply_cigar(const string &cigar, string &to_inflate,
            const char inflation_symbol) {
  std::istringstream iss(cigar);

  string inflated_seq;
  size_t n;
  char op;
  size_t i = 0;
  auto to_inflate_beg = std::begin(to_inflate);
  while (iss >> n >> op) {
    if (consumes_reference(op) && consumes_query(op)) {
      inflated_seq.append(to_inflate_beg + i, to_inflate_beg + i + n);
      i += n;
    }
    else if (consumes_query(op)) {
      // no addition of symbols to query
      i += n;
    }
    else if (consumes_reference(op)) {
      inflated_seq.append(n, inflation_symbol);
      // no increment of index within query
    }
  }

  // sum of total M/I/S/=/X/N operations must equal length of seq
  const size_t orig_len = to_inflate.length();
  if (i != orig_len)
    throw runtime_error("inconsistent number of qseq ops in cigar: " +
                        to_inflate + " "  + cigar + " " +
                        to_string(i) + " " +
                        to_string(orig_len));
  to_inflate.swap(inflated_seq);
}
