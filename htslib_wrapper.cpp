/*    Copyright (C) 2018 Andrew D. Smith
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
 */

#include "htslib_wrapper.hpp"

#include "GenomicRegion.hpp"
#include "MappedRead.hpp"

using std::string;

HTSFile::HTSFile(const std::string &fn) : is_good(true) {
  the_file = hts_open(fn.c_str(), "r");
  if (!the_file)
    throw std::runtime_error("cannot open file: " + fn);
  the_header = sam_hdr_read(the_file);
}

HTSFile::~HTSFile() {
  hts_close(the_file);
  bam_hdr_destroy(the_header);
}

static bool
get_next_read(HTSFile &h, bam1_t *b) {
  int rd = 0;
  // iterate to look for the next mapped read
  while ((rd = sam_read1(h.the_file, h.the_header, b)) >= 0 && b->core.tid < 0);
  return rd >= 0;
}

static void
bam1_to_GenomicRegion(const bam_hdr_t *h, const bam1_t *b,
                      GenomicRegion &r) {
  r.set_chrom(h->target_name[b->core.tid]);
  r.set_start(b->core.pos);
  r.set_end(bam_endpos(b));
  r.set_name(bam_get_qname(b));
  r.set_score(b->core.qual);
  r.set_strand(bam_is_rev(b) ? '-' : '+');
}

HTSFile &
operator>>(HTSFile &h, GenomicRegion &r) {
  bam1_t *b = bam_init1();
  if (get_next_read(h, b))
    bam1_to_GenomicRegion(h.the_header, b, r);
  else h.is_good = false;
  bam_destroy1(b);
  return h;
}


HTSFile &
operator>>(HTSFile &h, MappedRead &mr) {
  bam1_t *b = bam_init1();
  if (get_next_read(h, b)) {
    bam1_to_GenomicRegion(h.the_header, b, mr.r);
    const size_t readlen = b->core.l_qseq;
    mr.seq.resize(readlen);
    mr.scr.resize(readlen);
    unsigned char *qseq = bam_get_seq(b);
    unsigned char *qual = bam_get_qual(b);
    for (size_t i = 0; i < readlen; ++i) {
      // seq_nt16_str is declared in htslib/hts.h and decodes
      // nucleotides from nybble format
      mr.seq[i] = seq_nt16_str[bam_seqi(qseq, i)];
      mr.scr[i] = qual[i] + 33; // magic from htslib
    }
  }
  else h.is_good = false;
  bam_destroy1(b);
  return h;
}
