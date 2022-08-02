/* Part of Smith lab software
 *
 * Copyright (C) 2020 University of Southern California and
 *                    Andrew D. Smith
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


#ifndef BISULFITE_UTILS_HPP
#define BISULFITE_UTILS_HPP

#include <string>
#include <random>

#include "sam_record.hpp"

namespace bsflags {
  // ADS: this is our addition to the SAM flags, using the "free" bits
  /* 4096 0x1000 read is A-rich
   */
  static const uint16_t read_is_a_rich = 0x1000;
}

inline bool
is_t_rich(const sam_rec &sr) {
  return !samflags::check(sr.flags, bsflags::read_is_a_rich);
}

inline bool
is_a_rich(const sam_rec &sr) {
  return samflags::check(sr.flags, bsflags::read_is_a_rich);
}

inline void
set_t_rich(sam_rec &sr) {
  samflags::unset(sr.flags, bsflags::read_is_a_rich);
}

inline void
set_a_rich(sam_rec &sr) {
  samflags::set(sr.flags, bsflags::read_is_a_rich);
}

void
bisulfite_treatment(std::mt19937 &generator, std::string &seq,
                    double bs_rate = 1.0, double meth_rate = 0.0);

#endif
