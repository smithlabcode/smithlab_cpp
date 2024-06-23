/* Part of SMITHLAB software
 *
 * Copyright (C) 2008-2018 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include <algorithm>
#include <cmath>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_map>

#include "QualityScore.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>

using std::begin;
using std::ios_base;
using std::runtime_error;
using std::string;
using std::unordered_map;
using std::vector;

string strip_path(string full_path) {
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else
    ++start;
  return full_path.substr(start);
}

string strip_path_and_suffix(string full_path) {
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else
    ++start;
  size_t end = full_path.find_last_of('.');
  if (end == string::npos)
    end = full_path.length();
  return full_path.substr(start, end - start);
}

void parse_dir_basename_suffix(const string &full_path, string &dirname,
                               string &base_name, string &suffix) {
  const size_t base_index = full_path.find_last_of("/\\");
  if (base_index != string::npos)
    dirname = full_path.substr(0, base_index);
  else
    dirname = "";

  const size_t suffix_index = full_path.find_last_of(".");
  if (suffix_index != string::npos &&
      (base_index == string::npos || suffix_index > base_index))
    suffix = string(begin(full_path) + suffix_index + 1, end(full_path));
  else
    suffix = "";

  base_name = string(begin(full_path) + base_index + 1,
                     begin(full_path) + suffix_index);
}

void format_dir_basename_suffix(const string &dn, const string &bn,
                                const string &sf, string &full_path) {
  full_path = dn;
  if (!dn.empty() && bn.back() != '/')
    full_path += '/';

  full_path += bn;

  if (!sf.empty())
    full_path += (sf[0] == '.' || full_path.back() == '.') ? sf : "." + sf;
}

bool isdir(const char *filename) {
  struct stat buffer;
  stat(filename, &buffer);
  return S_ISDIR(buffer.st_mode);
}

bool is_fastq(const string &filename) {
  std::ifstream f(filename.c_str());
  char c = '\0';
  f >> c;
  f.close();
  return (c == '@');
}

////////////////////////////////////////////////////////////////////////
// Stuff dealing with FASTA format sequence files

bool is_valid_filename(const string &name, const string &filename_suffix) {
  const string suffix(name.substr(name.find_last_of(".") + 1));
  return (suffix == filename_suffix);
}

string path_join(const string &a, const string &b) {
  if (b.empty() || b[0] == '/')
    throw runtime_error("cannot prepend dir to file \"" + b + "\"");
  if (!a.empty() && a[a.length() - 1] == '/')
    return a + b;
  else
    return a + "/" + b;
}

void identify_chromosomes(const string &chrom_file, const string fasta_suffix,
                          unordered_map<string, string> &chrom_files) {
  vector<string> the_files;
  if (isdir(chrom_file.c_str())) {
    read_dir(chrom_file, fasta_suffix, the_files);
    for (size_t i = 0; i < the_files.size(); ++i)
      chrom_files[strip_path_and_suffix(the_files[i])] = the_files[i];
  }
  else
    chrom_files[strip_path_and_suffix(chrom_file)] = chrom_file;
}

void identify_and_read_chromosomes(const string &chrom_file,
                                   const string fasta_suffix,
                                   unordered_map<string, string> &chrom_files) {
  vector<string> the_files;
  if (isdir(chrom_file.c_str())) {
    read_dir(chrom_file, fasta_suffix, the_files);
  }
  else
    the_files.push_back(chrom_file);

  for (size_t i = 0; i < the_files.size(); ++i) {
    vector<string> names, seqs;
    read_fasta_file_short_names(the_files[i], names, seqs);
    for (size_t j = 0; j < names.size(); ++j)
      chrom_files[names[j]] = the_files[i];
  }
}

void read_dir(const string &dirname, string filename_suffix,
              vector<string> &filenames) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str())))
    throw runtime_error("could not open directory: " + dirname);

  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    if (is_valid_filename(ent->d_name, filename_suffix))
      filenames.push_back(path_join(dirname, string(ent->d_name)));
    errno = 0;
  }
  // check for some errors
  if (errno)
    throw runtime_error("error reading directory: " + dirname);
  if (filenames.empty())
    throw runtime_error("no valid files found in: " + dirname);
  closedir(dir);
}

bool is_sequence_line(const char *buffer) { return isvalid(buffer[0]); }

void parse_score_line(const char *buffer, vector<char> &scr) {
  for (const char *i = buffer; *i != '\0'; ++i)
    scr.push_back(*i);
}

inline bool is_fastq_name_line(size_t line_count) {
  return ((line_count & 3ul) == 0ul);
}

inline bool is_fastq_sequence_line(size_t line_count) {
  return ((line_count & 3ul) == 1ul);
}

inline bool is_fastq_score_name_line(size_t line_count) {
  return ((line_count & 3ul) == 2ul);
}

inline bool is_fastq_score_line(size_t line_count) {
  return ((line_count & 3ul) == 3ul);
}

void read_fastq_file(const char *filename, vector<string> &names,
                     vector<string> &sequences,
                     vector<vector<double>> &scores) {

  static const size_t INPUT_BUFFER_SIZE = 1000000;

  std::ifstream in(filename);
  if (!in)
    throw runtime_error("cannot open input file " + string(filename));

  string s, name;
  vector<char> scr;
  vector<vector<char>> scrs;
  bool first_line = true;
  // ADS: preprocessor stuff below is because is_sequence_line is only
  // used with asserts; consider removing variable
#ifndef NDEBUG
  bool is_sequence_line = false;
#endif
  size_t line_count = 0;
  while (!in.eof()) {
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw runtime_error("Line in " + name +
                          "\nexceeds max length: " + toa(INPUT_BUFFER_SIZE));
    if (in.gcount() == 0)
      break;

    // correct for dos carriage returns before newlines
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';

    if (is_fastq_name_line(line_count)) {
      if (buffer[0] != '@')
        throw runtime_error("invalid FASTQ name line: " + string(buffer));
      if (first_line == false && s.length() > 0) {
        names.push_back(name);
        sequences.push_back(s);
        scrs.push_back(scr);
      }
      else
        first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("@ "));
      s = "";
      scr.clear();
#ifndef NDEBUG
      is_sequence_line = true;
#endif
    }
    if (is_fastq_sequence_line(line_count)) {
      assert(is_sequence_line);
      s += buffer;
#ifndef NDEBUG
      is_sequence_line = false;
#endif
    }
    if (is_fastq_score_name_line(line_count)) {
      if (buffer[0] != '+')
        throw runtime_error("invalid FASTQ score name line: " + string(buffer));
    }
    if (is_fastq_score_line(line_count)) {
      parse_score_line(buffer, scr);
    }
    ++line_count;
  }
  if (!first_line && s.length() > 0) {
    names.push_back(name);
    sequences.push_back(s);
    scrs.push_back(scr);
  }

  bool phred_scores = true, solexa_scores = true;
  for (size_t i = 0; i < scrs.size() && phred_scores && solexa_scores; ++i) {
    phred_scores =
        phred_scores && (find_if(begin(scrs[i]), end(scrs[i]), [](char c) {
                           return !valid_phred_score(c);
                         }) == end(scrs[i]));
    solexa_scores =
        solexa_scores && (find_if(begin(scrs[i]), end(scrs[i]), [](char c) {
                            return !valid_solexa_score(c);
                          }) == end(scrs[i]));
  }

  if (!phred_scores && !solexa_scores)
    throw runtime_error("invalid quality scores in FASTQ file");

  for (size_t i = 0; i < scrs.size(); ++i) {
    scores.push_back(vector<double>(scrs[i].size()));
    for (size_t j = 0; j < scrs[i].size(); ++j)
      scores[i][j] = (solexa_scores)
                         ? quality_character_to_solexa(scrs[i][j] - 5)
                         : quality_character_to_phred(scrs[i][j]);
    scrs[i].clear();
  }
}

void read_fastq_file(const char *filename, vector<string> &names,
                     vector<string> &sequences, vector<string> &scores) {

  static const size_t INPUT_BUFFER_SIZE = 1000000;

  std::ifstream in(filename);
  if (!in)
    throw runtime_error("cannot open input file " + string(filename));

  string s, name, scr;
  bool first_line = true;
  // ADS: preprocessor stuff below is because is_sequence_line is only
  // used with asserts; consider removing variable
#ifndef NDEBUG
  bool is_sequence_line = false;
#endif
  size_t line_count = 0;
  while (!in.eof()) {
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw runtime_error("Line in " + name +
                          "\nexceeds max length: " + toa(INPUT_BUFFER_SIZE));
    if (in.gcount() == 0)
      break;

    // correct for dos carriage returns before newlines
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';

    if (is_fastq_name_line(line_count)) {
      if (buffer[0] != '@')
        throw runtime_error("invalid FASTQ name line: " + string(buffer));
      if (first_line == false && s.length() > 0) {
        names.push_back(name);
        sequences.push_back(s);
        scores.push_back(scr);
      }
      else
        first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("@ "));
#ifndef NDEBUG
      is_sequence_line = true;
#endif
    }
    if (is_fastq_sequence_line(line_count)) {
      assert(is_sequence_line);
      s = buffer;
#ifndef NDEBUG
      is_sequence_line = false;
#endif
    }
    if (is_fastq_score_name_line(line_count)) {
      if (buffer[0] != '+')
        throw runtime_error("invalid FASTQ score name line: " + string(buffer));
    }
    if (is_fastq_score_line(line_count)) {
      scr = buffer;
    }
    ++line_count;
  }
  if (!first_line && s.length() > 0) {
    names.push_back(name);
    sequences.push_back(s);
    scores.push_back(scr);
  }
}

void read_fasta_file_short_names(const string &filename, vector<string> &names,
                                 vector<string> &sequences) {

  std::ifstream in(filename);
  if (!in)
    throw runtime_error("cannot open input file " + filename);

  names.clear();
  sequences.clear();

  string line;
  while (getline(in, line)) {
    if (line[0] == '>') {
      const auto first_space = line.find_first_of(" \t", 1);
      if (first_space == string::npos)
        names.push_back(line.substr(1));
      else
        names.push_back(string(begin(line) + 1,
                               begin(line) + line.find_first_of(" \t", 1)));
      sequences.push_back(string());
    }
    else
      sequences.back() += line;
  }
}

void read_fasta_file(const string &filename, vector<string> &names,
                     vector<string> &sequences) {

  std::ifstream in(filename.c_str());
  if (!in)
    throw runtime_error("cannot open input file " + filename);

  names.clear();
  sequences.clear();

  string line;
  while (getline(in, line)) {

    if (line[0] == '>') {
      names.push_back(line.substr(1));
      sequences.push_back(string());
    }
    else
      sequences.back() += line;
  }
}

void read_fasta_file(const string &filename, const string &target,
                     string &sequence) {

  // read the sequence with the given name from a fasta file

  sequence = "";
  std::ifstream in(filename.c_str(), std::ios::binary);
  if (!in)
    throw runtime_error("cannot open input file " + filename);

  static const size_t INPUT_BUFFER_SIZE = 1000000;

  string s, name;

  bool first_line = true;
  while (!in.eof()) {
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw runtime_error("Line in " + name +
                          "\nexceeds max length: " + toa(INPUT_BUFFER_SIZE));
    // correct for dos carriage returns before newlines
    if (buffer[strlen(buffer) - 1] == '\r')
      buffer[strlen(buffer) - 1] = '\0';
    if (buffer[0] == '>') {
      if (first_line == false && s.length() > 0 && name == target) {
        std::swap(sequence, s);
        in.close();
        return;
      }
      else
        first_line = false;
      name = buffer;
      name = name.substr(name.find_first_not_of("> "));
      const size_t first_whitespace = name.find_first_of(" \t");
      if (first_whitespace != std::string::npos)
        name = name.substr(0, first_whitespace);
      s = "";
    }
    else if (name == target)
      s += buffer;
    in.peek();
  } // while

  if (!first_line && s.length() > 0 && name == target)
    std::swap(sequence, s);
}

void read_filename_file(const char *filename, vector<string> &filenames) {

  static const size_t INPUT_BUFFER_SIZE = 1000000;

  std::ifstream in(filename);
  if (!in)
    throw runtime_error("cannot open input file " + string(filename));

  while (!in.eof()) {
    char buffer[INPUT_BUFFER_SIZE + 1];
    in.getline(buffer, INPUT_BUFFER_SIZE);
    if (in.gcount() == static_cast<int>(INPUT_BUFFER_SIZE))
      throw runtime_error("Line in " + string(filename) +
                          "\nexceeds max length: " + toa(INPUT_BUFFER_SIZE));
    filenames.push_back(buffer);
    in.peek();
  }
}

size_t get_filesize(string filename) {
  std::ifstream f(filename.c_str());
  if (!f.good())
    return 0;

  size_t begin_pos = f.tellg();
  f.seekg(0, ios_base::end);
  size_t end_pos = f.tellg();
  f.close();

  return end_pos - begin_pos;
}

string basename(string filename) {
  const string s(filename.substr(0, filename.find_last_of(".")));
  const size_t final_slash = s.find_last_of("/");
  if (final_slash != string::npos)
    return s.substr(final_slash + 1);
  else
    return s;
}

void read_dir(const string &dirname, vector<string> &filenames) {
  DIR *dir;
  if (!(dir = opendir(dirname.c_str())))
    throw "could not open directory: " + dirname;

  errno = 0;
  struct dirent *ent;
  while ((ent = readdir(dir))) {
    filenames.push_back(path_join(dirname, string(ent->d_name)));
    errno = 0;
  }
  // check for some errors
  if (errno)
    throw "error reading directory: " + dirname;
  if (filenames.empty())
    throw "no valid files found in: " + dirname;
  closedir(dir);
}

bool is_valid_output_file(const string &filename) {
  // ADS: seems like there is no way around "access" and apparently
  // access is not a great solution anyway.
  if (std::filesystem::exists(filename))
    return (!std::filesystem::is_directory(filename) &&
            access(filename.c_str(), W_OK) == 0);
  else {
    // ADS: check if dir exists and is writeable
    // first get directory part
    string base_name, dir_part, suffix;
    parse_dir_basename_suffix(filename, dir_part, base_name, suffix);
    if (dir_part.empty())
      dir_part = "./";
    // check if directory part exists and is writeable
    return (access(dir_part.c_str(), F_OK | W_OK) == 0);
  }
}

bool has_gz_ext(const string &filename) {
  const string ext(".gz");
  return filename.size() >= ext.size() &&
         filename.compare(filename.size() - ext.size(), ext.size(), ext) == 0;
}
