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

#include <iterator>

#include "SAM.hpp"
#include "smithlab_utils.hpp"
#include "MappedRead.hpp"

#include <fstream>

using std::string;
using std::vector;

SAM::SAM(string mapper) : mapper(mapper) {}

SAM::SAM(string mapper, string line) : mapper(mapper) {
  std::istringstream iss(line);
  if (mapper.compare("bsmap") == 0) {
    if (!(iss >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
        >> mate_name >> mate_start >> seg_len >> seq >> qual
        >> mismatch_str >> strand_str))
      throw SMITHLABException("malformed line in bsmap SAM format:\n" + line);

    IS_TRICH = flag & 0x40;
    get_mismatch_bsmap();
  }
  else if (mapper.compare("bismark") == 0) {
    if (!(iss >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
        >> mate_name >> mate_start >> seg_len >> seq >> qual
        >> edit_distance_str >> mismatch_str >> meth_call_str
        >> read_conv_str >> genome_conv_str))
      throw SMITHLABException("malformed line in bismark SAM format:\n" + line);

    IS_TRICH = read_conv_str.substr(read_conv_str.size() - 2)
      == genome_conv_str.substr(genome_conv_str.size() - 2);
    get_mismatch_bismark();
  }
  else if (mapper.compare("bs_seeker") == 0) {
    if (!(iss >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
        >> mate_name >> mate_start >> seg_len >> seq >> qual))
      throw SMITHLABException("malformed line in bs_seeker SAM format:\n" + line);

    // seems bs_seeker doesn't keep mismatch information.
    mismatch = 0;

    // bs_seeker also doesn't keep sequencing quality information?
    string new_qual(seq.size(), 'h');
    qual = new_qual;
  }
  else if (mapper.compare("unknown") == 0) {
    if (!(iss >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
        >> mate_name >> mate_start >> seg_len >> seq >> qual
        ))
      throw SMITHLABException("malformed line in SAM format:\n" + line);

    mismatch = 0;
  }
  else {
    throw SMITHLABException("Unsupported mapper for SAM format: " + mapper
        + ". Please try using 'unknown'.");
  }
}

MappedRead
SAM::GetMappedRead()  {
  MappedRead mr;
  //SAM is 1-based
  GenomicRegion r(chrom,start - 1,start - 1 + seq.size());
  
  if (mapper.compare("bismark") == 0)
    get_mr_bismark(mr,r);
  else if (mapper.compare("bsmap") == 0)
    get_mr_bsmap(mr,r);
  else if (mapper.compare("bs_seeker") == 0)
    get_mr_bsseeker(mr,r);
  else
    get_mr_general(mr,r);

  return mr;
}

void
SAM::get_mr_bsmap(MappedRead &mr, GenomicRegion &r) const {
  // bsmap stores the original seq of each read regardless of which strand
  // it mapped to.
  string new_seq, new_qual;
  apply_CIGAR(new_seq, new_qual);

  if (is_revcomp())
    r.set_strand('-');
  else
    r.set_strand('+');

  r.set_score(mismatch);
  mr.r = r;
  mr.r.set_end(r.get_start() + new_seq.size()); //update region length
  mr.seq = new_seq;
  mr.scr = new_qual;

  string strand, bs_forward;
  get_strand(strand, bs_forward);

  if (is_pairend()) {
    if (is_Trich()) {
      mr.r.set_name(name + "/1");
    }
    else {
      mr.r.set_name(name + "/2");
    }
    assert((is_Trich() && bs_forward == "+")
      || (is_Arich() && bs_forward == "-"));
  }
  else {
    //single end
    assert(bs_forward == "+");
    assert(is_revcomp() == (strand != bs_forward));
    mr.r.set_name(name);
  }
  mr.r.set_strand(strand[0]);
}

void
SAM::get_mr_bismark(MappedRead &mr, GenomicRegion &r) {
  
  // string new_seq, new_qual;
  // apply_CIGAR(new_seq, new_qual);

  // if (is_revcomp())
  //   r.set_strand('-');
  // else
  //   r.set_strand('+');

  // // if a read is mapped to - strand, bismark stores the + strand seq of
  // // reference genome rather than the seq of that read. I'm not sure
  // // about the orientation of CIGAR string in bismark. But we do need the
  // // original sequence in .mr
  // if (is_revcomp())
  // {
  //   revcomp_inplace(new_seq);
  //   std::reverse(new_qual.begin(), new_qual.end());
  // }

  // if (CIGAR.find_first_of("IDSHPX") != string::npos)
  // {
  //   throw SMITHLABException("Only support bowtie 1 reaults without gaps");
  // }
  
  const string read_conv_mode = read_conv_str.substr(read_conv_str.size() - 2);
  const string genome_conv_mode =
      genome_conv_str.substr(genome_conv_str.size() - 2);
  if (genome_conv_mode == "CT")
  {
    r.set_strand('+');
    
    string new_seq, new_qual;
    if (read_conv_mode == "CT")
      apply_CIGAR(new_seq, new_qual);
    else
    {
      revcomp_inplace(seq);
      std::reverse(qual.begin(), qual.end());
      apply_CIGAR(new_seq, new_qual);
      revcomp_inplace(new_seq);
      std::reverse(new_qual.begin(), new_qual.end());
      revcomp_inplace(seq);
      std::reverse(qual.begin(), qual.end());
    }
    
    r.set_score(mismatch);
    mr.r = r;
    mr.r.set_end(r.get_start() + new_seq.size()); //update region length
    mr.seq = new_seq;
    mr.scr = new_qual;
    mr.r.set_name(name);

//    std::cout << mr << std::endl;
    
  }
  else if (genome_conv_mode == "GA") 
  {
    r.set_strand('-');

    string new_seq, new_qual;
    if (read_conv_mode == "GA")
      apply_CIGAR(new_seq, new_qual);
    else
    {
      revcomp_inplace(seq);
      std::reverse(qual.begin(), qual.end());
      apply_CIGAR(new_seq, new_qual);
      revcomp_inplace(new_seq);
      std::reverse(new_qual.begin(), new_qual.end());
      revcomp_inplace(seq);
      std::reverse(qual.begin(), qual.end());
    }

    std::transform(new_seq.begin(), new_seq.end(), new_seq.begin(), complement);
    r.set_score(mismatch);
    mr.r = r;
    mr.r.set_end(r.get_start() + new_seq.size()); //update region length
    mr.seq = new_seq;
    mr.scr = new_qual;
    mr.r.set_name(name);

//    std::cout << mr << std::endl;

  }
}

void
SAM::get_mr_bsseeker(MappedRead &mr, GenomicRegion &r) const {
  throw SMITHLABException("NOT FULLY SUPPORTED");

  // string new_seq, new_qual;

  // apply_CIGAR(new_seq, new_qual);

  // if (is_revcomp())
  //   r.set_strand('-');
  // else
  //   r.set_strand('+');

  // // if a read is mapped to - strand, bs_seeker stores the + strand seq of
  // // reference genome rather than the seq of that read. I'm not sure
  // // about the orientation of CIGAR string in bs_seeker. But we do need the
  // // original sequence in .mr
  // if (is_revcomp())
  // {
  //   revcomp_inplace(new_seq);
  //   std::reverse(new_qual.begin(), new_qual.end());
  // }

  // r.set_score(mismatch);
  // mr.r = r;
  // mr.r.set_end(r.get_start() + new_seq.size()); //update region length
  // mr.seq = new_seq;
  // mr.scr = new_qual;

  // mr.r.set_name(name);
}

void
SAM::get_mr_general(MappedRead &mr, GenomicRegion &r) const {
  string new_seq, new_qual;
  apply_CIGAR(new_seq, new_qual);

  if (is_revcomp())
    r.set_strand('-');
  else
    r.set_strand('+');

  r.set_score(mismatch);
  mr.r = r;
  mr.r.set_end(r.get_start() + new_seq.size()); //update region length
  mr.seq = new_seq;
  mr.scr = new_qual;

  mr.r.set_name(name);
}

bool
SAM::load_read_from_line(std::istream& the_stream) {
  if (mapper.compare("bsmap") == 0) {
    if (!(the_stream >> name >> flag >> chrom >> start
          >> mapq_score >> CIGAR >> mate_name >> mate_start
          >> seg_len >> seq >> qual >> mismatch_str >> strand_str))
      return false;

    get_mismatch_bsmap();
    return true;
  }
  else if (mapper.compare("bismark") == 0) {
    if (!(the_stream >> name >> flag >> chrom >> start
          >> mapq_score >> CIGAR >> mate_name
          >> mate_start >> seg_len >> seq >> qual
          >> edit_distance_str >> mismatch_str >> meth_call_str
          >> read_conv_str >> genome_conv_str))
      return false;

    get_mismatch_bismark();
    return true;
  }
  else if (mapper.compare("bs_seeker") == 0) {
    if (!(the_stream >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
        >> mate_name >> mate_start >> seg_len >> seq >> qual))
      return false;

    // seems bs_seeker doesn't keep mismatch information.
    mismatch = 0;

    // bs_seeker also doesn't keep sequencing quality information?
    string new_qual(seq.size(), 'h');
    qual = new_qual;
    return true;
  }
  else if (mapper.compare("unknown") == 0) {
    if (!(the_stream >> name >> flag >> chrom >> start >> mapq_score >> CIGAR
        >> mate_name >> mate_start >> seg_len >> seq >> qual))
      return false;

    mismatch = 0;
    return true;
  }
  
  return false;
}

void
SAM::apply_CIGAR(string &new_seq, string &new_qual) const {
    assert(seq.size() == qual.size());
    assert(new_seq.size() == 0 && new_qual.size() == 0);
    size_t n;
    char op;
    size_t i = 0;

    std::istringstream iss(CIGAR);
    while (iss >> n >> op)
    {
        switch (op)
        {
        case 'M':
            new_seq += seq.substr(i, n);
            new_qual += qual.substr(i, n);
            i += n;
            break;
        case 'I':
            i += n;
            break;
        case 'D':
            new_seq += string(n, 'N');
            new_qual += string(n, 'B');
            break;
        case 'S':
            i += n;
            break;
        case 'H':
            ;
            break;
        case 'P':
            ;
            break;
        case '=':
            ;
            break;
        case 'X':
            ;
            break;
        }
    }
    
    assert(i == seq.length());
    assert(new_seq.size() == new_qual.size());
}

void
SAM::get_mismatch_bsmap() {
  mismatch = atoi(mismatch_str.substr(5).c_str());
}

void
SAM::get_mismatch_bismark() {
  /*
  the result of this function might not be accurate, because if a sequencing
  error occurs on a cytosine, then it probably will be reported as a convertion
  */
  size_t edit_distance;
  edit_distance = atoi(edit_distance_str.substr(5).c_str());

  int convert_count = 0;
  const char *temp = meth_call_str.substr(5).c_str();
  while(*temp != '\0') {
    if (*temp == 'x' || *temp == 'h' || *temp == 'z')
      ++convert_count;
    ++temp;
  }

  mismatch = edit_distance - convert_count;
}

void
SAM::get_strand(string &strand, string &bs_forward) const {
  strand = strand_str.substr(5, 1);
  bs_forward = strand_str.substr(6, 1);
  if (bs_forward == "-") strand = strand == "+" ? "-" : "+";
}

std::istream& 
operator>>(std::istream& the_stream, SAM &r) {
  //To-do: do not load unmapped reads?
  if (!r.load_read_from_line(the_stream))
      //throw SMITHLABException("Unable to load read from input file.");
      the_stream.setstate(std::ios::badbit);
    
    char c;
    while ((c = the_stream.get()) != '\n' && the_stream);
    
    if (c != '\n')
      //throw SMITHLABException("Unable to load read from input file.");
      the_stream.setstate(std::ios::badbit);
    
    // the_stream.peek();
    if (the_stream.eof())
      //throw SMITHLABException("Unable to load read from input file.");
      the_stream.setstate(std::ios::badbit);

  return the_stream;
}

std::ostream& 
operator<<(std::ostream& the_stream, SAM &r) {
  //output as .mr format
  MappedRead mr = r.GetMappedRead();

  return the_stream << mr;
}
