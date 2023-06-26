/* Part of Smith lab software
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

#ifndef SAM_RECORD_HPP
#define SAM_RECORD_HPP

#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iterator>
#include <cstdint>

// from 30 April 2020 SAM documentation
// 1    0x1   template having multiple segments in sequencing
// 2    0x2   each segment properly aligned according to the aligner
// 4    0x4   segment unmapped
// 8    0x8   next segment in the template unmapped
// 16   0x10  SEQ being reverse complemented
// 32   0x20  SEQ of the next segment in the template being reverse complemented
// 64   0x40  the first segment in the template
// 128  0x80  the last segment in the template
// 256  0x100 secondary alignment
// 512  0x200 not passing filters, such as platform/vendor quality controls
// 1024 0x400 PCR or optical duplicate
// 2048 0x800 supplementary alignment

namespace samflags {

  // ADS: names of flags adjusted to how we typically interpret
  static const uint16_t read_paired = 0x1;
  static const uint16_t read_pair_mapped = 0x2;
  static const uint16_t read_unmapped = 0x4;
  static const uint16_t mate_unmapped = 0x8;
  static const uint16_t read_rc = 0x10;
  static const uint16_t mate_rc = 0x20;
  static const uint16_t template_first = 0x40;
  static const uint16_t template_last = 0x80;
  static const uint16_t secondary_aln = 0x100;
  static const uint16_t below_quality = 0x200;
  static const uint16_t pcr_duplicate = 0x400;
  static const uint16_t supplementary_aln = 0x800;

  constexpr bool
  check(const uint16_t to_check, const uint16_t &f) {return to_check & f;}
  inline void
  set(uint16_t &to_set, const uint16_t f) {to_set |= f;}
  inline void
  unset(uint16_t &to_unset, const uint16_t f) {to_unset &= ~f;}
}

std::string
format_sam_flags(const uint16_t the_flags);

// Col Field Type Regexp/Range Brief description
// 1 QNAME Query template NAME
// 2 FLAG bitwise FLAG
// 3 RNAME Reference sequence NAME11
// 4 POS 1-based leftmost mapping POSition
// 5 MAPQ MAPping Quality
// 6 CIGAR CIGAR string
// 7 RNEXT Reference name of the mate/next read
// 8 PNEXT Position of the mate/next read
// 9 TLEN observed Template LENgth
//10 SEQ segment SEQuence
//11 QUAL Phred-scaled base QUALity+33

class sam_rec {
public:
  // ADS: instance vars *are* in SAM order
  std::string qname;
  uint16_t flags;
  std::string rname;
  uint32_t pos;
  uint8_t mapq;
  std::string cigar;
  std::string rnext;
  uint32_t pnext;
  int32_t tlen;
  std::string seq;
  std::string qual;
  std::vector<std::string> tags;
  sam_rec() : flags(0), pos(0), mapq(255), pnext(0), tlen(0) {}
  explicit sam_rec(const std::string &line);
  sam_rec(const std::string &_qname,
          const uint16_t _flags,
          const std::string &_rname,
          const uint32_t &_pos,
          const uint8_t &_mapq,
          const std::string &_cigar,
          const std::string &_rnext,
          const uint32_t &_pnext,
          const int32_t &_tlen,
          const std::string &_seq,
          const std::string &_qual) :
    qname(_qname),
    flags(_flags),
    rname(_rname),
    pos(_pos),
    mapq(_mapq),
    cigar(_cigar),
    rnext(_rnext),
    pnext(_pnext),
    tlen(_tlen),
    seq(_seq),
    qual(_qual) {}
  void add_tag(const std::string &the_tag) {tags.push_back(the_tag);}
  bool empty() const;
  size_t estimate_line_size() const;
  std::string tostring() const;
};

inline bool
check_flag(const sam_rec &sr, const uint16_t f) {
  return samflags::check(sr.flags, f);
}

inline void
set_flag(sam_rec &sr, const uint16_t f) {
  return samflags::set(sr.flags, f);
}

inline void
unset_flag(sam_rec &sr, const uint16_t f) {
  return samflags::unset(sr.flags, f);
}

std::istream &
operator>>(std::istream &in, sam_rec &r);

std::ostream &
operator<<(std::ostream &the_stream, const sam_rec &r);

void
inflate_with_cigar(const sam_rec &sr, std::string &to_inflate,
                   const char inflation_symbol = 'N');

// assumes program name is at argv[0]
template<typename T>
void
write_pg_line(const int argc, const char **argv,
              const std::string &program_name,
              const std::string &program_version, T &out) {
  out << "@PG\t";

  // empty program, should never happen
  if (argc == 0) {
    out << '\n';
  }
  else {
    if (!program_name.empty())
      out << "ID:" << program_name << '\t';
    if (!program_version.empty())
      out << "VN:" << program_version << '\t';

    std::ostringstream the_command;
    copy(argv, argv + argc,
         std::ostream_iterator<const char*>(the_command, " "));
    out << "CL:\"" << the_command.str() << "\"" << std::endl;
  }
}

template<typename T>
static void
write_sam_header(const std::vector<std::string> &chrom_names,
                 const std::vector<T> &chrom_sizes,
                 const std::string program_name,
                 const std::string program_version,
                 const int argc, const char **argv,
                 std::ostream &out) {
  static const std::string SAM_VERSION = "1.0";

  // sam version
  out <<"@HD" << '\t' << "VN:" << SAM_VERSION << '\n'; // sam version

  // chromosome sizes
  const size_t n_chroms = chrom_names.size();
  for (size_t i = 0; i < n_chroms; ++i) {
    out << "@SQ" << '\t'
        << "SN:" << chrom_names[i] << '\t'
        << "LN:" << chrom_sizes[i] << '\n';
  }

  // program details
  out << "@PG" << '\t'
      << "ID:" << program_name << '\t'
      << "VN:" << program_version << '\t';

  // how the program was run
  std::ostringstream the_command;
  copy(argv, argv + argc, std::ostream_iterator<const char*>(the_command, " "));
  out << "CL:\"" << the_command.str() << "\"" << std::endl;
}

#endif
