/*
 *    Part of SMITHLAB software
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

#include "bisulfite_utils.hpp"
#include <cstdlib>
#include <iostream>
#include <random>

using std::string;

void bisulfite_treatment(std::mt19937 &generator, string &seq, double bs_rate,
                         double meth_rate) {

  std::uniform_real_distribution<double> unif(0.0, 1.0);

  const size_t seq_len = seq.length() - 1;
  for (size_t i = 0; i < seq_len; ++i)
    if (toupper(seq[i]) == 'C') {
      // CpG
      if (toupper(seq[i + 1]) == 'G' &&
          (unif(generator) > meth_rate || unif(generator) < bs_rate)) {
        seq[i] = 'T';
      }
      // Regular C
      else if (unif(generator) < bs_rate)
        seq[i] = 'T';
    }
  if (toupper(seq[seq_len]) == 'C' && unif(generator) < bs_rate)
    seq[seq_len] = 'T';
}
