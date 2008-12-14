/*
 *    Part of RMAP software
 *
 *    Copyright (C) 2008 Cold Spring Harbor Laboratory, 
 *                       University of Southern California and
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

#include "FastReadQuality.hpp"

#include <cmath>
#include <iomanip>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

size_t FastReadQuality::score_mask = 0;
size_t FastReadQuality::segments = 0;
size_t FastReadQuality::read_width = 0;
size_t FastReadQuality::right_most_bit = 0;

double FastReadQuality::scaler = 0;
double FastReadQuality::min_quality_score = 0;
double FastReadQuality::max_quality_score = 0;

////////////////////////////////////////////////////////////////////////
// WORD PAIR

double
FastReadQuality::value_to_quality(size_t val) {
  return ((double(val)/scaler)*(max_quality_score - min_quality_score) + 
	  read_width*min_quality_score)/max_quality_score;
}

size_t
FastReadQuality::quality_to_value(double quality) {
  return static_cast<size_t>(scaler*std::min(1.0, (quality - min_quality_score)/(max_quality_score - 
										 min_quality_score)));
}

size_t
FastReadQuality::Words::quality_to_value(double quality) {
  return static_cast<size_t>(scaler*std::min(1.0, (quality - min_quality_score)/(max_quality_score - min_quality_score)));
}

double
FastReadQuality::Words::value_to_quality(size_t val) {
  return min_quality_score + (val/scaler)*(max_quality_score - min_quality_score);
}

FastReadQuality::Words::Words(const vector<vector<double> > &s) : 
  a_vec(0), c_vec(0), g_vec(0), t_vec(0) {
  const vector<vector<double> >::const_iterator limit = s.end();
  for (vector<vector<double> >::const_iterator i(s.begin()); i != limit; ++i) {
    a_vec = ((a_vec << n_val_bits) + (quality_to_value((*i)[0])));
    c_vec = ((c_vec << n_val_bits) + (quality_to_value((*i)[1])));
    g_vec = ((g_vec << n_val_bits) + (quality_to_value((*i)[2])));
    t_vec = ((t_vec << n_val_bits) + (quality_to_value((*i)[3])));
  }
  if (s.size()*n_val_bits < rmap_bits::word_size) {
    const size_t additional_shift = (rmap_bits::word_size - s.size()*n_val_bits);
    a_vec <<= additional_shift;
    c_vec <<= additional_shift;
    g_vec <<= additional_shift;
    t_vec <<= additional_shift;
  }
}


FastReadQuality::Words::Words(string::const_iterator i, const string::const_iterator limit) :
  a_vec(0), c_vec(0), g_vec(0), t_vec(0) {
  const size_t length = limit - i;
  for (; i < limit; ++i) {
    a_vec = ((a_vec << n_val_bits) | contains_a(base2int(*i)));
    c_vec = ((c_vec << n_val_bits) | contains_c(base2int(*i)));
    g_vec = ((g_vec << n_val_bits) | contains_g(base2int(*i)));
    t_vec = ((t_vec << n_val_bits) | contains_t(base2int(*i)));
  }
  if (length*n_val_bits < rmap_bits::word_size) {
    const size_t additional_shift = (rmap_bits::word_size - length*n_val_bits);
    a_vec <<= additional_shift;
    c_vec <<= additional_shift;
    g_vec <<= additional_shift;
    t_vec <<= additional_shift;
  }
}


size_t
FastReadQuality::Words::get_val(MASK_t mask, MASK_t base_vec, size_t pos) {
  // 00 -> A, 01 -> C, 10 -> G, 11 -> T
  const MASK_t selector = (FastReadQuality::low_val_bits << (pos - 1)*n_val_bits);
  return (((mask & base_vec) & selector) >> (pos - 1)*n_val_bits);
}


string
FastReadQuality::Words::bits2string(size_t mask, size_t bits) {
  string s;
  size_t selector = rmap_bits::high_bit;
  for (size_t i = 0; i < rmap_bits::word_size; ++i) {
    s += (selector & bits & mask) ? '1' : '0';
    selector >>= 1;
  }
  return s;
}


string
FastReadQuality::Words::tostring_bits(size_t mask) const {
  return (bits2string(mask, a_vec) + "\n" +
	  bits2string(mask, c_vec) + "\n" +
	  bits2string(mask, g_vec) + "\n" +
	  bits2string(mask, t_vec) + "\n");

}


string
FastReadQuality::Words::tostring_values(size_t mask) const {
  std::ostringstream ss;
  for (size_t i = FastReadQuality::segment_size; i > 0; --i)
    ss << std::setw(3) << get_val(mask, a_vec, i) << " ";
  ss << endl;
  for (size_t i = FastReadQuality::segment_size; i > 0; --i)
    ss << std::setw(3) << get_val(mask, c_vec, i) << " ";
  ss << endl;
  for (size_t i = FastReadQuality::segment_size; i > 0; --i)
    ss << std::setw(3) << get_val(mask, g_vec, i) << " ";
  ss << endl;
  for (size_t i = FastReadQuality::segment_size; i > 0; --i)
    ss << std::setw(3) << get_val(mask, t_vec, i) << " ";
  ss << endl;

  return ss.str();
}


////////////////////////////////////////////////////////////////////////
// FAST READ

void
FastReadQuality::set_read_properties(const size_t rw, const double minqs, const double maxqs) {
  read_width = rw;
  segments = static_cast<size_t>(std::ceil(rw*n_val_bits/static_cast<float>(rmap_bits::word_size)));
  right_most_bit = (rmap_bits::word_size - (rw*n_val_bits % rmap_bits::word_size));
  score_mask = (rmap_bits::all_ones << right_most_bit);
  
  min_quality_score = minqs;
  max_quality_score = maxqs;
  scaler = pow(2.0, n_val_bits) - 1;
}

FastReadQuality::FastReadQuality(const vector<vector<double> > &s) {
  for (size_t i = 0; i < segments - 1; ++i) {
    const vector<vector<double> > this_seg(s.begin() + i*segment_size,
					   s.begin() + (i + 1)*segment_size);
    words.push_back(Words(this_seg));
  }
  const vector<vector<double> > this_seg(s.begin() + (segments - 1)*segment_size, s.end());
  words.push_back(Words(this_seg));
}

FastReadQuality::FastReadQuality(const string &s) {
  for (size_t i = 0; i < segments - 1; ++i)
    words.push_back(Words(s.begin() + i*segment_size, s.begin() + (i + 1)*segment_size));
  words.push_back(Words(s.begin() + (segments - 1)*segment_size, s.end()));
}

string
FastReadQuality::tostring_values() const {
  std::ostringstream ss;
  for (size_t i = 0; i < words.size() - 1; ++i)
    ss << words[i].tostring_values(rmap_bits::all_ones) << endl;
  ss << words.back().tostring_values(score_mask);
  return ss.str();
}

string
FastReadQuality::tostring_bits() const {
  std::ostringstream ss;
  for (size_t i = 0; i < words.size() - 1; ++i)
    ss << words[i].tostring_bits(rmap_bits::all_ones) << endl;
  ss << words.back().tostring_bits(score_mask);
  return ss.str();
}