/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2011 University of Southern California and
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

#include "GenomicRegion.hpp"

#include <exception>

#include <cassert>
#include <fstream>
#include <unordered_map>

using std::ostringstream;
using std::runtime_error;
using std::string;
using std::unordered_map;
using std::vector;

unordered_map<string, chrom_id_type> SimpleGenomicRegion::fw_table_in;
unordered_map<chrom_id_type, string> SimpleGenomicRegion::fw_table_out;

chrom_id_type SimpleGenomicRegion::assign_chrom(const std::string &c) {
  auto chr_id(fw_table_in.find(c));
  if (chr_id == fw_table_in.end()) {
    const chrom_id_type r = fw_table_in.size();
    fw_table_in[c] = r;
    fw_table_out[r] = c;
    // return r;
  }
  // else return chr_id->second;
  return fw_table_in[c];
}

string SimpleGenomicRegion::retrieve_chrom(chrom_id_type i) {
  auto chr_name(fw_table_out.find(i));
  // assert(chr_name != fw_table_out.end());
  return chr_name->second;
}

SimpleGenomicRegion::SimpleGenomicRegion(const GenomicRegion &r)
    : chrom(assign_chrom(r.get_chrom())), start(r.get_start()),
      end(r.get_end()) {}

SimpleGenomicRegion::SimpleGenomicRegion(const char *s, const size_t len) {
  size_t i = 0;

  // the chrom
  while (isspace(s[i]) && i < len)
    ++i;
  size_t j = i;
  while (!isspace(s[i]) && i < len)
    ++i;
  chrom = assign_chrom(string(s + j, i - j));

  // start of the region (a positive integer)
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  while (!isspace(s[i]) && i < len)
    ++i;
  start = atoi(s + j);

  // end of the region (a positive integer)
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  while (!isspace(s[i]) && i < len)
    ++i;
  end = atoi(s + j);
}

string SimpleGenomicRegion::tostring() const {
  std::ostringstream s;
  s << retrieve_chrom(chrom) << "\t" << start << "\t" << end;
  return s.str();
}

bool SimpleGenomicRegion::contains(const SimpleGenomicRegion &other) const {
  return chrom == other.chrom && start <= other.start && other.end <= end;
}

bool SimpleGenomicRegion::overlaps(const SimpleGenomicRegion &other) const {
  return chrom == other.chrom &&
         ((start < other.end && other.end <= end) ||
          (start <= other.start && other.start < end) || other.contains(*this));
}

size_t SimpleGenomicRegion::distance(const SimpleGenomicRegion &other) const {
  if (chrom != other.chrom)
    return std::numeric_limits<size_t>::max();
  else if (overlaps(other) || other.overlaps(*this))
    return 0;
  else
    return (end < other.start) ? other.start - end + 1 : start - other.end + 1;
}

bool SimpleGenomicRegion::operator<(const SimpleGenomicRegion &rhs) const {
  return (get_chrom() < rhs.get_chrom() ||
          (chrom == rhs.chrom &&
           (start < rhs.start || (start == rhs.start && (end < rhs.end)))));
}

bool SimpleGenomicRegion::less1(const SimpleGenomicRegion &rhs) const {
  return (get_chrom() < rhs.get_chrom() ||
          (chrom == rhs.chrom &&
           (end < rhs.end || (end == rhs.end && start < rhs.start))));
}

bool SimpleGenomicRegion::operator<=(const SimpleGenomicRegion &rhs) const {
  return !(rhs < *this);
}

bool SimpleGenomicRegion::operator==(const SimpleGenomicRegion &rhs) const {
  return (chrom == rhs.chrom && start == rhs.start && end == rhs.end);
}

bool SimpleGenomicRegion::operator!=(const SimpleGenomicRegion &rhs) const {
  return (chrom != rhs.chrom || start != rhs.start || end != rhs.end);
}

#include <iostream>

unordered_map<string, chrom_id_type> GenomicRegion::fw_table_in;
unordered_map<chrom_id_type, string> GenomicRegion::fw_table_out;

chrom_id_type GenomicRegion::assign_chrom(const std::string &c) {
  auto chr_id(fw_table_in.find(c));
  if (chr_id == fw_table_in.end()) {
    const chrom_id_type r = fw_table_in.size();
    fw_table_in[c] = r;
    fw_table_out[r] = c;
    return r;
  }
  else
    return chr_id->second;
}

string GenomicRegion::retrieve_chrom(chrom_id_type i) {
  unordered_map<chrom_id_type, string>::const_iterator chr_name(
      fw_table_out.find(i));
  assert(chr_name != fw_table_out.end());
  return chr_name->second;
}

GenomicRegion::GenomicRegion(const char *s, const size_t len) {
  size_t i = 0;

  // the chrom
  while (isspace(s[i]) && i < len)
    ++i;
  if (i == len)
    throw runtime_error(
        "malformatted BED file contains only one"
        "column in the line below "
        "(a properly formatted BED file must contain at least three):\n" +
        string(s));

  size_t j = i;
  while (!isspace(s[i]) && i < len)
    ++i;
  if (i == len)
    throw runtime_error(
        "malformatted BED file contains only two "
        "columns in the line below "
        "(a properly formatted BED file must contain at least three):\n" +
        string(s));

  chrom = assign_chrom(string(s + j, i - j));

  // start of the region (a positive integer)
  while (isspace(s[i]) && i < len)
    ++i;

  j = i;
  start = atoi(s + j);
  while (!isspace(s[i]) && i < len)
    ++i;

  // end of the region (a positive integer)
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  end = atoi(s + j);
  while (!isspace(s[i]) && i < len)
    ++i;

  // name of the region
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  while (!isspace(s[i]) && i < len)
    ++i;
  name = string(s + j, i - j);

  // score of the region (floating point)
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  score = atof(s + j);
  while (!isspace(s[i]) && i < len)
    ++i;

  // strand
  while (isspace(s[i]) && i < len)
    ++i;
  j = i;
  strand = *(s + j);
  while (!isspace(s[i]) && i < len)
    ++i;

  // ADS: This is a hack!!!
  if (strand != '-')
    strand = '+';
}

string GenomicRegion::tostring() const {
  std::ostringstream s;
  s << retrieve_chrom(chrom) << "\t" << start << "\t" << end;
  if (!name.empty())
    s << "\t" << name << "\t" << score << "\t" << strand;
  return s.str();
}

bool GenomicRegion::contains(const GenomicRegion &other) const {
  return chrom == other.chrom && start <= other.start && other.end <= end;
}

bool GenomicRegion::overlaps(const GenomicRegion &other) const {
  return chrom == other.chrom &&
         ((start < other.end && other.end <= end) ||
          (start <= other.start && other.start < end) || other.contains(*this));
}

size_t GenomicRegion::distance(const GenomicRegion &other) const {
  if (chrom != other.chrom)
    return std::numeric_limits<size_t>::max();
  else if (overlaps(other) || other.overlaps(*this))
    return 0;
  else
    return (end < other.start) ? other.start - end + 1 : start - other.end + 1;
}

bool GenomicRegion::operator<(const GenomicRegion &rhs) const {
  return ((chrom == rhs.chrom &&
           (start < rhs.start ||
            (start == rhs.start &&
             (end < rhs.end || (end == rhs.end &&
                                (strand < rhs.strand
                                 // || (strand == rhs.strand && name < rhs.name)
                                 )))))) ||
          get_chrom() < rhs.get_chrom());
}

bool GenomicRegion::less1(const GenomicRegion &rhs) const {
  return ((chrom == rhs.chrom &&
           (end < rhs.end || (end == rhs.end &&
                              (start < rhs.start ||
                               (start == rhs.start &&
                                (strand < rhs.strand
                                 // || (strand == rhs.strand && name < rhs.name)
                                 )))))) ||
          get_chrom() < rhs.get_chrom());
}

bool GenomicRegion::operator<=(const GenomicRegion &rhs) const {
  return !(rhs < *this);
}

bool GenomicRegion::operator==(const GenomicRegion &rhs) const {
  return (chrom == rhs.chrom && start == rhs.start && end == rhs.end &&
          name == rhs.name && score == rhs.score && strand == rhs.strand);
}

bool GenomicRegion::operator!=(const GenomicRegion &rhs) const {
  return (chrom != rhs.chrom || start != rhs.start || end != rhs.end ||
          name != rhs.name || score != rhs.score || strand != rhs.strand);
}

void separate_chromosomes(
    const vector<SimpleGenomicRegion> &regions,
    vector<vector<SimpleGenomicRegion>> &separated_by_chrom) {
  typedef unordered_map<chrom_id_type, vector<SimpleGenomicRegion>> Separator;
  Separator separator;
  for (vector<SimpleGenomicRegion>::const_iterator i = regions.begin();
       i != regions.end(); ++i) {
    const chrom_id_type the_chrom(i->chrom);
    if (separator.find(the_chrom) == separator.end())
      separator[the_chrom] = vector<SimpleGenomicRegion>();
    separator[the_chrom].push_back(*i);
  }
  separated_by_chrom.clear();
  for (Separator::iterator i = separator.begin(); i != separator.end(); ++i)
    separated_by_chrom.push_back(i->second);
}

void separate_chromosomes(const vector<GenomicRegion> &regions,
                          vector<vector<GenomicRegion>> &separated_by_chrom) {
  typedef unordered_map<chrom_id_type, vector<GenomicRegion>> Separator;
  Separator separator;
  for (vector<GenomicRegion>::const_iterator i = regions.begin();
       i != regions.end(); ++i) {
    const chrom_id_type the_chrom(i->chrom);
    if (separator.find(the_chrom) == separator.end())
      separator[the_chrom] = vector<GenomicRegion>();
    separator[the_chrom].push_back(*i);
  }
  separated_by_chrom.clear();
  for (Separator::iterator i = separator.begin(); i != separator.end(); ++i)
    separated_by_chrom.push_back(i->second);
}

static inline auto is_header_line(const string &line) -> bool {
  static constexpr auto browser_label = "browser";
  static constexpr auto browser_label_len = 7u;
  for (auto i = 0u; i < browser_label_len; ++i)
    if (line[i] != browser_label[i])
      return false;
  return true;
}

static inline auto is_track_line(const string &line) -> bool {
  static constexpr auto track_label = "track";
  static constexpr auto track_label_len = 5u;
  for (auto i = 0u; i < track_label_len; ++i)
    if (line[i] != track_label[i])
      return false;
  return true;
}

void ReadBEDFile(const string &filename, vector<GenomicRegion> &the_regions) {
  std::ifstream in(filename);
  if (!in)
    throw runtime_error("failed to open file " + filename);

  string line;
  while (getline(in, line))
    if (!is_header_line(line) && !is_track_line(line))
      the_regions.push_back(GenomicRegion(line));
}

void ReadBEDFile(const string &filename,
                 vector<SimpleGenomicRegion> &the_regions) {
  std::ifstream in(filename);
  if (isdir(filename.c_str()))
    throw runtime_error("BED file is a directory: " + filename);

  if (!in.good())
    throw runtime_error("cannot open input file " + filename);

  string line;
  while (getline(in, line))
    if (!is_header_line(line) && !is_track_line(line))
      the_regions.push_back(SimpleGenomicRegion(line));
}
