/* Part of Smith Lab software
 *
 * Copyright (C) 2020 University of Southern California and
 *                    Guilherme De Sena Brandine and Andrew D. Smith
 *
 * Authors: Guilherme De Sena Brandine and Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include "sam_record.hpp"

#include <cstdint>
#include <regex>
#include <sstream>

#include "cigar_utils.hpp"
#include "smithlab_utils.hpp"

using std::begin;
using std::end;
using std::istream;
using std::istringstream;
using std::ostream;
using std::ostringstream;
using std::regex;
using std::runtime_error;
using std::string;
using std::to_string;

// ADS: this is for debugging purposes
string format_sam_flags(const uint16_t the_flags) {
  ostringstream oss;
  using samflags::check;
  oss << "read_paired: " << check(the_flags, samflags::read_paired) << '\n'
      << "read_pair_mapped: " << check(the_flags, samflags::read_pair_mapped)
      << '\n'
      << "read_unmapped: " << check(the_flags, samflags::read_unmapped) << '\n'
      << "mate_unmapped: " << check(the_flags, samflags::mate_unmapped) << '\n'
      << "read_rc: " << check(the_flags, samflags::read_rc) << '\n'
      << "mate_rc: " << check(the_flags, samflags::mate_rc) << '\n'
      << "template_first: " << check(the_flags, samflags::template_first)
      << '\n'
      << "template_last: " << check(the_flags, samflags::template_last) << '\n'
      << "secondary_aln: " << check(the_flags, samflags::secondary_aln) << '\n'
      << "below_quality: " << check(the_flags, samflags::below_quality) << '\n'
      << "pcr_duplicate: " << check(the_flags, samflags::pcr_duplicate) << '\n'
      << "supplementary_aln: " << check(the_flags, samflags::supplementary_aln);
  return oss.str();
}

// static bool
// valid_cigar(const string &cigar, const string &seq) {
//   return cigar_qseq_ops(cigar) == seq.size();
// }

// static bool
// valid_seq(const string &read) {
//   return regex_match(read, regex("\\*|[A-Za-z=.]+"));
// }

// static bool
// valid_qual(const string &qual) {
//   return regex_match(qual, regex("[!-~]+"));
// }

size_t sam_rec::estimate_line_size() const {
  static const size_t all_field_estimates = 100;
  return qname.size() + rname.size() + qual.size() + all_field_estimates;
}

string sam_rec::tostring() const {
  string out;
  out.reserve(estimate_line_size());
  out.append(qname + "\t" + to_string(flags) + "\t" + rname + "\t" +
             to_string(pos) + "\t" + to_string(static_cast<unsigned>(mapq)) +
             "\t" + cigar + "\t" + rnext + "\t" + to_string(pnext) + "\t" +
             to_string(tlen) + "\t" + seq + "\t" + qual);
  for (auto it(begin(tags)); it != end(tags); ++it)
    out.append("\t" + *it);

  return out;
}

ostream &operator<<(std::ostream &the_stream, const sam_rec &r) {
  the_stream << r.tostring();
  return the_stream;
}

sam_rec::sam_rec(const string &line) {
  /*
    istringstream iss; // ADS: change to set the buffer from "line"
    iss.rdbuf()->pubsetbuf(const_cast<char*>(line.c_str()), line.size());*/
  istringstream iss(line);      // ADS: unfortunate macos stuff?
  int32_t will_become_mapq = 0; // to not read mapq as character
                                // since it's uint8_t
  if (!(iss >> qname >> flags >> rname >> pos >> will_become_mapq >> cigar >>
        rnext >> pnext >> tlen >> seq >> qual))
    throw runtime_error("incorrect SAM record:\n" + line);
  if (will_become_mapq < 0 || will_become_mapq > 255)
    throw runtime_error("invalid mapq in SAM record: " + line);
  mapq = static_cast<uint8_t>(will_become_mapq);

  string tmp;
  while (iss >> tmp)
    tags.push_back(tmp);
}

bool sam_rec::empty() const { return pos == 0; }

void inflate_with_cigar(const sam_rec &sr, string &to_inflate,
                        const char inflation_symbol) {
  apply_cigar(sr.cigar, to_inflate, inflation_symbol);
}
