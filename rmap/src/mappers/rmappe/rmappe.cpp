/*    rmappe: a program for mapping paired-end Solexa reads
 *    Copyright (C) 2009 University of Southern California and
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

#include <string>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <cmath>

#include "FastRead.hpp"
#include "FastReadQuality.hpp"
#include "rmap_os.hpp"
#include "SeedMaker.hpp"
#include "MapResultPE.hpp"
#include "OptionParser.hpp"

#include <tr1/unordered_map>

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::ostream_iterator;
using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::make_pair;

using std::tr1::unordered_multimap;

typedef unordered_multimap<size_t, size_t> SeedHash;

void
get_read_matches(const size_t the_seed, const vector<string> &reads,
		 SeedHash &seed_hash) {
  for (size_t i = 0; i < reads.size(); ++i) {
    const size_t truncate_to = min(reads[i].length(), SeedMaker::max_seed_part);
    const string s(reads[i].substr(0, truncate_to));
    if (SeedMaker::valid_seed(the_seed, s)) {
      const size_t the_key = (the_seed & SeedMaker::make_read_word(s));
      seed_hash.insert(SeedHash::value_type(the_key, i));
    }
  }
}


template <class T> void
map_reads(const string &chrom, const size_t chrom_id, const size_t profile, 
	  const size_t read_width, const size_t max_diffs, 
	  const size_t min_sep, const size_t max_sep, 
	  const vector<T> &reads_left, const vector<T> &reads_right, 
	  const SeedHash &seed_hash, const bool strand, 
	  vector<MultiMapResultPE> &best_maps) {
  
  MASK_t bad_bases = rmap_bits::all_ones;
  MASK_t read_word = rmap_bits::all_zeros;
  T fast_read;
  
  const size_t key_diff = read_width - min(read_width, SeedMaker::max_seed_part);
  const size_t chrom_size = chrom.size();

  vector<T> possible_lefts(max_sep);
  
  // const size_t lookup_mask = table_size - 1;
  
  size_t chrom_offset = 0;
  for (chrom_offset = 0; chrom_offset < key_diff; ++chrom_offset)
    fast_read.shift(base2int(chrom[chrom_offset]));

  for (; chrom_offset < read_width - 1; ++chrom_offset) {
    const size_t key_base = base2int(chrom[chrom_offset - key_diff]);
    fast_read.shift(base2int(chrom[chrom_offset]));
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
  }
  
  for (; chrom_offset < chrom_size; ++chrom_offset) {
    const size_t base = base2int(chrom[chrom_offset]);
    const size_t key_base = base2int(chrom[chrom_offset - key_diff]);
    
    fast_read.shift(base);
    SeedMaker::update_bad_bases(key_base, bad_bases);
    SeedMaker::update_read_word(key_base, read_word);
    
    if ((bad_bases & profile) == 0) {
      std::pair<SeedHash::const_iterator, SeedHash::const_iterator>
	bucket(seed_hash.equal_range((read_word & profile)));
      if (bucket.first != bucket.second) {
	
	const SeedHash::const_iterator  limit(bucket.second);
	for (SeedHash::const_iterator to_test(bucket.first); to_test != limit; ++to_test) {

	  const size_t score = reads_right[to_test->second].score(fast_read);

	  if (score <= max_diffs) {

	    const size_t lookback_limit = chrom_offset - min_sep;
	    for (size_t i = chrom_offset - max_sep; i < lookback_limit; ++i) {
	      
	      const size_t pair_score = score + 
		reads_left[to_test->second].score(possible_lefts[i % max_sep]);
	      
	      const vector<MultiMapResultPE>::iterator current(best_maps.begin() + 
							       to_test->second);
	      if (pair_score <= current->score) {
		const size_t left_start = i - read_width + 1;
		current->add(score, chrom_id, left_start, 
			     chrom_offset - read_width + 1, strand);
	      }
	    }
	  }
	}
      }
      possible_lefts[chrom_offset % max_sep] = fast_read;
    }
  }
}


bool
good_read(const vector<vector<double> > &read) {
  size_t bad_count = 0;
  for (size_t i = 0; i < read.size(); ++i)
    bad_count += (*std::max_element(read[i].begin(), read[i].end()) <= 0);
  return (bad_count <= 2);
}


void
clean_reads(const size_t min_match_score, const size_t left_width,
	    const vector<string> &input_read_names, 
	    vector<string> &reads_left, vector<string> &reads_right, 
	    vector<vector<vector<double> > > &scores_left,
	    vector<vector<vector<double> > > &scores_right,
	    vector<vector<string> > &read_names) {
  vector<vector<vector<double> > > good_scores_left, good_scores_right;
  vector<string> good_reads_left, good_reads_right;
  for (size_t i = 0; i < reads_left.size(); ++i) {
    if (good_read(scores_left[i]) && good_read(scores_right[i])) {
      read_names.push_back(vector<string>(1, input_read_names[i]));
      good_reads_left.push_back(reads_left[i]);
      good_reads_right.push_back(reads_right[i]);
      good_scores_left.push_back(scores_right[i]);
      good_scores_right.push_back(scores_right[i]);
    }
  }
  reads_left.swap(good_reads_left);
  reads_right.swap(good_reads_right);
  scores_left.swap(good_scores_left);
  scores_right.swap(good_scores_right);
}


void
sites_to_regions(const size_t read_width, const vector<string> &chrom, 
		 const vector<size_t> &chrom_sizes, 
		 const vector<vector<string> > &read_names, 
		 vector<MultiMapResultPE> &bests, const size_t max_diffs, 
		 vector<GenomicRegion> &hits) {
  for (size_t i = 0; i < bests.size(); ++i)
    if (!bests[i].mr.empty()) {
      sort(bests[i].mr.begin(), bests[i].mr.end());
      for (size_t j = 0; j < bests[i].mr.size(); ++j)
	if (j == 0 || bests[i].mr[j - 1] < bests[i].mr[j]) {
	  const size_t chrom_id = bests[i].mr[j].chrom;
	  const size_t start = bests[i].mr[j].strand ? 
	    bests[i].mr[j].site : 
	    chrom_sizes[chrom_id] - bests[i].mr[j].site2 - read_width;
	  const size_t end = bests[i].mr[j].strand ? 
	    bests[i].mr[j].site2 + read_width : 
	    chrom_sizes[chrom_id] - bests[i].mr[j].site;
	  const size_t score = bests[i].mr[j].score;
	  const char strand = ((bests[i].mr[j].strand) ? '+' : '-');
	  for (size_t k = 0; k < read_names[i].size(); ++k)
	    hits.push_back(GenomicRegion(chrom[chrom_id], start, end,
					 read_names[i][k], score, strand));
	}
    }
}


static char
to_base_symbol(char c) {return (isvalid(c)) ? toupper(c) : 'N';}


void
clean_reads(const size_t max_diffs, const size_t read_width,
	    const vector<string> &input_read_names,
	    vector<string> &reads_left, vector<string> &reads_right, 
	    vector<vector<string> > &read_names) {
  vector<pair<pair<string, string>, string> > sorter;
  for (size_t i = 0; i < reads_left.size(); ++i) {
    transform(reads_left[i].begin(), reads_left[i].end(), reads_left[i].begin(),
	      std::ptr_fun(&to_base_symbol));
    transform(reads_right[i].begin(), reads_right[i].end(), reads_right[i].begin(),
	      std::ptr_fun(&to_base_symbol));
    if (count(reads_left[i].begin(), reads_left[i].end(), 'N') <=
	static_cast<int>(max_diffs) &&
	count(reads_right[i].begin(), reads_right[i].end(), 'N') <=
	static_cast<int>(max_diffs))
      sorter.push_back(make_pair(make_pair(reads_left[i], reads_right[i]),
				 input_read_names[i]));
  }
  sort(sorter.begin(), sorter.end());
  reads_left.clear();
  reads_right.clear();
  for (size_t i = 0; i < sorter.size(); ++i) {
    if (i == 0 || sorter[i - 1].first != sorter[i].first) {
      reads_left.push_back(sorter[i].first.first);
      reads_right.push_back(sorter[i].first.second);
      read_names.push_back(vector<string>(1, sorter[i].second));
    }
    else read_names.back().push_back(sorter[i].second);
  }
}


void
sites_to_regions(const size_t read_width, const vector<string> &chrom, 
		 const vector<size_t> &chrom_sizes, 
		 const vector<vector<string> > &read_names, 
		 vector<MultiMapResultPE> &bests, const size_t min_match_score, 
		 const double max_quality_score,
		 vector<GenomicRegion> &hits) {
  for (size_t i = 0; i < bests.size(); ++i)
    if (!bests[i].mr.empty()) {
      sort(bests[i].mr.begin(), bests[i].mr.end());
      for (size_t j = 0; j < bests[i].mr.size(); ++j)
	if (j == 0 || bests[i].mr[j - 1] < bests[i].mr[j]) {
	  const size_t chrom_id = bests[i].mr[j].chrom;
	  const size_t start = bests[i].mr[j].strand ? 
	    bests[i].mr[j].site : 
	    chrom_sizes[chrom_id] - bests[i].mr[j].site - read_width;
	  const size_t end = start + read_width;
	  const size_t score = bests[i].mr[j].score;
	  const char strand = ((bests[i].mr[j].strand) ? '+' : '-');
	  for (size_t k = 0; k < read_names[i].size(); ++k)
	    hits.push_back(GenomicRegion(chrom[chrom_id], start, end,
					 read_names[i][k], 
					 FastReadQuality::value_to_quality(score), 
					 strand));
	}
    }
}


void
write_non_uniques(string filename, const vector<pair<string, size_t> > &ambigs) {
  std::ofstream out(filename.c_str());
  for (size_t i = 0; i < ambigs.size(); ++i)
    out << ambigs[i].first << "\t" << ambigs[i].second << endl;
  out.close();
}


template <class T> void
eliminate_ambigs(const size_t max_mismatches,
		 vector<MultiMapResultPE> &best_maps, vector<vector<string> > &read_names,
		 vector<string> &reads_left, vector<string> &reads_right, 
		 vector<pair<string, size_t> > &ambigs,
		 vector<T> &fast_reads_left, vector<T> &fast_reads_right) {
  
  size_t j = 0;
  for (size_t i = 0; i < best_maps.size(); ++i) {
    best_maps[i].collapse();
    if (best_maps[i].ambiguous() && best_maps[i].score <= max_mismatches)
      for (size_t k = 0; k < read_names[i].size(); ++k)
	ambigs.push_back(make_pair(read_names[i][k], best_maps[i].score));
    else {
      best_maps[j] = best_maps[i];
      read_names[j].swap(read_names[i]);
      reads_left[j].swap(reads_left[i]);
      reads_right[j].swap(reads_right[i]);
      fast_reads_left[j] = fast_reads_left[i];
      fast_reads_right[j] = fast_reads_right[i];
      ++j;
    }
  }
  best_maps.erase(best_maps.begin() + j, best_maps.end());
  read_names.erase(read_names.begin() + j, read_names.end());
  reads_left.erase(reads_left.begin() + j, reads_left.end());
  reads_right.erase(reads_right.begin() + j, reads_right.end());
  fast_reads_left.erase(fast_reads_left.begin() + j, fast_reads_left.end());
  fast_reads_right.erase(fast_reads_right.begin() + j, fast_reads_right.end());
}


void
process_score_data(const size_t read_width, const size_t max_mismatches,
		   double &max_quality_score, double &max_match_score,
		   vector<vector<vector<double> > > &scores_left,
		   vector<vector<vector<double> > > &scores_right) {
  double min_quality_score = numeric_limits<double>::max();
  for (size_t i = 0; i < scores_left.size(); ++i)
    for (size_t j = 0; j < scores_left[i].size(); ++j) {
      max_quality_score = 
	max(max_quality_score, 
	    max(*max_element(scores_left[i][j].begin(), scores_left[i][j].end()),
		*max_element(scores_right[i][j].begin(), scores_right[i][j].end())));
      min_quality_score = 
	min(min_quality_score, 
	    min(*min_element(scores_left[i][j].begin(), scores_left[i][j].end()),
		*min_element(scores_right[i][j].begin(), scores_right[i][j].end())));
    }
  
  FastReadQuality::set_read_properties(read_width, 0, max_quality_score - 
				       min_quality_score);
  const double max_final_score = 
    max_mismatches*(max_quality_score - min_quality_score);
  max_match_score = max_mismatches*
    FastReadQuality::quality_to_value(max_final_score/max_mismatches);
  
  for (size_t i = 0; i < scores_left.size(); ++i) {
    for (size_t j = 0; j < scores_left[i].size(); ++j) {
      const double max_score_left = *max_element(scores_left[i][j].begin(), 
						 scores_left[i][j].end());
      const double max_score_right = *max_element(scores_right[i][j].begin(), 
						  scores_right[i][j].end());
      for (size_t k = 0; k < scores_left[i][j].size(); ++k) {
	scores_left[i][j][k] = max_score_left - scores_left[i][j][k];
	scores_right[i][j][k] = max_score_right - scores_right[i][j][k];
      }
    }
  }
}


void
split_scores(const size_t left_width,
	     vector<vector<vector<double> > > &scores, 
	     vector<vector<vector<double> > > &left,
	     vector<vector<vector<double> > > &right) {
  for (size_t i = 0; i < scores.size(); ++i) {
    left.push_back(vector<vector<double> >(scores[i].begin(), 
					   scores[i].begin() + left_width));
    right.push_back(vector<vector<double> >(scores[i].begin() + left_width,
					    scores[i].begin() + 2*left_width));
    scores[i].clear();
  }
  scores.clear();
}


void
split_reads(const size_t left_width, vector<string> &reads,
	    vector<string> &left, vector<string> &right) {
  for (size_t i = 0; i < reads.size(); ++i) {
    left.push_back(reads[i].substr(0, left_width));
    right.push_back(reads[i].substr(left_width, left_width));
  }
  reads.clear();
}


int 
main(int argc, const char **argv) {
  try {
    
    string chrom_file;
    string filenames_file;
    string outfile;
    string prb_file;
    string ambiguous_file;
    string fasta_suffix = "fa";
    
    size_t n_seeds = 0;
    size_t seed_weight = 0;
    size_t read_width = 0;
    size_t max_mismatches = 0;
    size_t max_mappings = 1;

    size_t left_width = 0;
    size_t max_sep = 200;
    size_t min_sep = 0;
    
    bool VERBOSE = false;
    
    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse("rmappe", "The rmap mapping tool for "
			   "Solexa reads (paired-end version)");
    opt_parse.add_opt("output", 'o', "Name of output file (default: stdout)", 
		      false , outfile);
    opt_parse.add_opt("chrom", 'c', "FASTA file or dir containing chromosome(s)", 
		      true , chrom_file);
    opt_parse.add_opt("suffix", 's', "suffix of FASTA files (assumes -c indicates dir)", 
		      false , fasta_suffix);
    opt_parse.add_opt("filenames", 'F', "file listing names of chromosome files", 
		      false , filenames_file);
    opt_parse.add_opt("prb", 'p', "file with quality scores (prb format)", 
		      false, prb_file);
    opt_parse.add_opt("seeds", 'S', "number of seeds", true , n_seeds);
    opt_parse.add_opt("hit", 'h', "weight of hit", true , seed_weight);
    opt_parse.add_opt("width", 'w', "width of reads", false, read_width);
    opt_parse.add_opt("mismatch", 'm', "maximum allowed mismatches", 
		      true , max_mismatches);
    opt_parse.add_opt("ambiguous", 'a', "file to write names of ambiguously mapped reads", 
		      false , ambiguous_file);
    opt_parse.add_opt("min-sep", '\0', "min separation between ends", 
		      false, min_sep);
    opt_parse.add_opt("max-sep", '\0', "max separation between ends", 
		      false, max_sep);
    opt_parse.add_opt("max-map", 'M', "maximum allowed mappings for a read", 
		      false, max_mappings);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string reads_file = leftover_args.front();
    /****************** END COMMAND LINE OPTIONS *****************/

    const bool USING_QUALITY = !prb_file.empty();
      
    // Get the reads
    vector<string> input_read_names, reads;
    read_fasta_file(reads_file.c_str(), input_read_names, reads);
    
    vector<vector<vector<double> > > scores;
    if (USING_QUALITY) {
      read_prb_file(prb_file.c_str(), scores);
      if (scores.size() != reads.size())
	throw RMAPException("different number of reads in prb and reads files");
    }
    
    if (VERBOSE)
      cerr << "TOTAL READS:    " << reads.size() << endl;

    if (read_width == 0) {
      read_width = reads.front().size();
      for (size_t i = 1; i < reads.size(); ++i)
	read_width = min(read_width, reads[i].length());
    }
    
    if (VERBOSE)
      cerr << "READ WIDTH:     " << read_width << endl;
    
    if (left_width > read_width/2) {
      throw RMAPException("incorrect read width");
    }
    else if (left_width == 0)
      left_width = std::floor(read_width/2);
    
    vector<string> reads_left, reads_right;
    split_reads(left_width, reads, reads_left, reads_right);
    vector<vector<vector<double> > > scores_left, scores_right;
    if (USING_QUALITY)
      split_scores(left_width, scores, scores_left, scores_right);
    
    // prepare the reads
    vector<vector<string> > read_names;
    if (USING_QUALITY)
      clean_reads(max_mismatches, left_width, input_read_names, 
		  reads_left, reads_right, scores_left, scores_right, 
		  read_names);
    else clean_reads(max_mismatches, left_width, input_read_names,
		     reads_left, reads_right, read_names);
    
    if (VERBOSE)
      cerr << "READS AFTER QC: " << reads_left.size() << endl;
    
    // convert the reads into FastReads
    vector<FastRead> fast_reads_left, fast_reads_right;
    vector<FastReadQuality> fast_reads_left_q, fast_reads_right_q;
    double max_match_score = 0, max_quality_score = 0;
    if (USING_QUALITY) {
      process_score_data(read_width, max_mismatches,
			 max_quality_score, max_match_score, 
			 scores_left, scores_right);
      if (VERBOSE)
	cerr << "MAX_MATCH_SCORE=" << max_match_score << endl
	     << "MAX_QUALITY_SCORE=" << max_quality_score << endl;
      for (size_t i = 0; i < scores_left.size(); ++i) {
	fast_reads_left_q.push_back(FastReadQuality(scores_left[i]));
	fast_reads_right_q.push_back(FastReadQuality(scores_right[i]));
	scores_left[i].clear();
	scores_right[i].clear();
      }
      scores_left.clear();
      scores_right.clear();
    }
    else {
      FastRead::set_read_width(left_width);
      for (size_t i = 0; i < reads_left.size(); ++i) {
	fast_reads_left.push_back(FastRead(reads_left[i]));
	fast_reads_right.push_back(FastRead(reads_right[i]));
      }
      // Should be able to do this: reads.clear();
    }
    
    vector<string> chrom_files;
    if (!filenames_file.empty())
      read_filename_file(filenames_file.c_str(), chrom_files);
    else if (isdir(chrom_file.c_str())) 
      read_dir(chrom_file, fasta_suffix, chrom_files);
    else chrom_files.push_back(chrom_file);
    
    if (VERBOSE) {
      cerr << endl << "chromosome files found (approx size):" << endl;
      for (vector<string>::const_iterator i = chrom_files.begin();
	   i != chrom_files.end(); ++i)
	cerr << *i << " (" << roundf(get_filesize(*i)/1e06) << "Mbp)" << endl;
    }
    
    vector<string> chrom_names;
    for (size_t i = 0; i < chrom_files.size(); ++i)
      chrom_names.push_back(basename(chrom_files[i]));
    
    if (VERBOSE)
      cerr << endl << "scanning chromosomes:" << endl;

    vector<size_t> the_seeds;
    SeedMaker::first_last_seeds(min(left_width, SeedMaker::max_seed_part),
				n_seeds, seed_weight, the_seeds);
    
    if (VERBOSE) {
      cerr << "SEED STRUCTURES" << endl;
      for (size_t i = 0; i < the_seeds.size(); ++i)
	cerr << bits2string_masked(rmap_bits::all_ones, the_seeds[i]) << endl;
      cerr << endl;
    }
    vector<size_t> chrom_sizes(chrom_files.size());
    vector<pair<string, size_t> > ambigs;

    MultiMapResultPE::init(max_mappings);
    vector<MultiMapResultPE> best_maps(reads_left.size(), MultiMapResultPE(max_mismatches));
    
    for (size_t j = 0; j < the_seeds.size() && !reads_left.empty(); ++j) {
      for (size_t i = 0; i < chrom_files.size() && !reads_left.empty(); ++i) {
	if (VERBOSE)
	  cerr << "[SEED:" << j + 1 << "/" << the_seeds.size() << "] [SEEDING] ";
	
	SeedHash seed_hash;
	get_read_matches(the_seeds[j], reads_right, seed_hash);
	
	vector<string> chrom_names, chrom;
	if (VERBOSE)
	  cerr << "[LOADING CHROM] ";
	read_fasta_file(chrom_files[i].c_str(), chrom_names, chrom);
	if (j == 0)
	  chrom_sizes[i] = chrom.front().length();
	
	if (VERBOSE)
	  cerr << "[SCANNING=" << chrom_names.front() << "] ";

	if (USING_QUALITY)
	  map_reads(chrom.front(), i, the_seeds[j], left_width, max_match_score,
		    min_sep, max_sep,
		    fast_reads_left_q, fast_reads_right_q,
		    seed_hash, true, best_maps);
	else map_reads(chrom.front(), i, the_seeds[j], left_width, max_mismatches,
		       min_sep, max_sep,
		       fast_reads_left, fast_reads_right,
		       seed_hash, true, best_maps);
	
	revcomp_inplace(chrom.front());

	if (USING_QUALITY)
	  map_reads(chrom.front(), i, the_seeds[j], left_width, max_match_score,
		    min_sep, max_sep,
		    fast_reads_left_q, fast_reads_right_q,
		    seed_hash, false, best_maps);
	else map_reads(chrom.front(), i, the_seeds[j], left_width, max_mismatches,
		       min_sep, max_sep,
		       fast_reads_left, fast_reads_right,
		       seed_hash, false, best_maps);
	
	if (VERBOSE)
	  cerr << "[CLEANING=" << chrom_names.front() << "] ";
	if (USING_QUALITY)
	  eliminate_ambigs(0, best_maps, read_names, reads_left, reads_right, 
			   ambigs, fast_reads_left_q, fast_reads_right_q);
	else
	  eliminate_ambigs(0, best_maps, read_names, reads_left, reads_right,
			   ambigs, fast_reads_left, fast_reads_right);
	if (VERBOSE)
	  cerr << "[AMBIG=" << ambigs.size() << "] " << endl;
      }
    }

    if (VERBOSE)
      cerr << endl << "[eliminating ambiguous reads...";
    
    if (USING_QUALITY)
      eliminate_ambigs(max_mismatches, best_maps, read_names, 
		       reads_left, reads_right, ambigs, 
		       fast_reads_left_q, fast_reads_right_q);
    else eliminate_ambigs(max_mismatches, best_maps, read_names, 
			  reads_left, reads_right, ambigs, 
			  fast_reads_left, fast_reads_right);
    
    if (!ambiguous_file.empty())
      write_non_uniques(ambiguous_file, ambigs);
    if (VERBOSE)
      cerr << "done]" << endl;
    
    // transform best matches into BED format
    vector<GenomicRegion> hits;

    if (USING_QUALITY)
      sites_to_regions(left_width, chrom_names, chrom_sizes,
		       read_names, best_maps, max_match_score, 
		       max_quality_score, hits);
    else sites_to_regions(left_width, chrom_names, chrom_sizes,
			  read_names, best_maps, max_mismatches, hits);
    
    // Output the results
    std::ostream* out = (!outfile.empty()) ? 
      new std::ofstream(outfile.c_str()) : &cout;
    copy(hits.begin(), hits.end(), ostream_iterator<GenomicRegion>(*out, "\n"));
    if (out != &cout) delete out;

  }
  catch (const RMAPException &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}