/*
 *    Part of SMITHLAB_CPP software
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Meng Zhou, Qiang Song, Andrew Smith
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

#ifndef SAM_HPP
#define SAM_HPP

#include "smithlab_utils.hpp"
#include "MappedRead.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <tr1/unordered_map>
#include <limits>

// All structs are written w/o consideration of header lines.
class SAM {
public:
  //constructors
  SAM() {}
  SAM(std::string mapper);
  SAM(std::string mapper, std::string line);

  //accessors
  std::string get_name() const {return name;}

  //flag status checkers
  bool is_pairend() const {return flag & 0x1;}
  bool is_singlend() const {return !(is_pairend());}
  bool is_mapping_paired() const {return flag & 0x2;}
  bool is_unmapped() const {return flag & 0x4;}
  bool is_mapped() const {return !(is_unmapped());}
  bool is_revcomp() const {return flag & 0x10;}
  bool is_Trich() const {return flag & 0x40;}
  bool is_Arich() const {return flag & 0x80;}
  bool is_secondary() const {return flag & 0x100;}
  bool is_primary() const {return !(is_secondary());}

  MappedRead GetMappedRead() const;
  bool load_read_from_line(std::istream& the_stream);

private:
  //standard SAM format fields
  std::string mapper;
  std::string name;
  std::string chrom;
  std::string CIGAR;
  std::string mate_name;
  std::string seq;
  std::string qual;

  //bismark specific fields
  std::string edit_distance_str;
  std::string mismatch_str;
  std::string meth_call_str;
  std::string read_conv_str;
  std::string genome_conv_str;

  //bsmap specific fields
  std::string strand_str;

  //some other (standard) attributes
  size_t flag;
  size_t start;
  size_t mapq_score;
  size_t mate_start;
  size_t mismatch;
  int seg_len;

  //functions to get MappedRead
  void get_mr_bsmap(MappedRead &mr, GenomicRegion &r) const;
  void get_mr_bismark(MappedRead &mr, GenomicRegion &r) const;
  void get_mr_bsseeker(MappedRead &mr, GenomicRegion &r) const;
  void get_mr_general(MappedRead &mr, GenomicRegion &r) const;

  //other functions associated with information extraction
  void apply_CIGAR(std::string &new_seq, std::string &new_qual) const;
  void get_mismatch_bsmap();
  void get_mismatch_bismark();
  void get_strand(std::string &strand, std::string &bs_forward) const;
};

std::istream& 
operator>>(std::istream& the_stream, SAM &r);

std::ostream& 
operator<<(std::ostream& the_stream, SAM &r);

#endif
