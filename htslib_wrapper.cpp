/* Part of SMITHLAB_CPP software
 *
 * Copyright (C) 2013-2019 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Meng Zhou, Qiang Song, Andrew Smith
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

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "htslib_wrapper.hpp"
#include "smithlab_utils.hpp"
#include "MappedRead.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::endl;
using std::runtime_error;

char check_htslib_wrapper() {return 1;}

SAMReader::SAMReader(const string &fn) :
  filename(fn), good(true) {

  if (!(hts = hts_open(filename.c_str(), "r")))
    throw runtime_error("cannot open file: " + filename);

  if (hts_get_format(hts)->category != sequence_data)
    throw runtime_error("file format appears wrong: " + filename);

  if (!(hdr = sam_hdr_read(hts)))
    throw runtime_error("failed to read header from file: " + filename);

  if (!(b = bam_init1()))
    throw runtime_error("failed to read record from file: " + filename);
}

SAMReader::~SAMReader() {
  if (hdr) {
    bam_hdr_destroy(hdr);
    hdr = 0;
  }
  if (b) {
    bam_destroy1(b);
    b = 0;
  }
  if (hts) {
    assert(hts_close(hts) >= 0);
    hts = 0;
  }
  good = false;
}


SAMReader&
operator>>(SAMReader &reader, sam_rec &aln) {
  reader.get_sam_record(aln);
  return reader;
}

/////////////////////////////////////////////
//// general facility for SAM format
/////////////////////////////////////////////

bool
SAMReader::get_sam_record(sam_rec &sr) {
  int rd_ret = 0;
  if ((rd_ret = sam_read1(hts, hdr, b)) >= 0) {
    int fmt_ret = 0;
    if ((fmt_ret = sam_format1(hdr, b, &hts->line)) <= 0)
      throw runtime_error("failed reading record from: " + filename);
    sr = sam_rec(hts->line.s);
    good = true; //reader.get_sam_record(reader.hts->line.s, sr);
    // ADS: line below seems to be needed when the file format is SAM
    // because the hts_getline for parsing SAM format files within
    // sam_read1 only get called when "(fp->line.l == 0)". For BAM
    // format, it does not seem to matter.
    hts->line.l = 0;
  }
  else if (rd_ret == -1)
    good = false;
  else // rd_ret < -1
    throw runtime_error("failed to read record from file: " + filename);
  return good;
}

string
SAMReader::get_header() const {
  return hdr->text; // includes newline
}
