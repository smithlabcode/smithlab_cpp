/*  Copyright (C) 2020 University of Southern California
 *                     and Andrew D. Smith
 *
 *  Authors: Andrew D. Smith
 *
 *  This is free software: you can redistribute it and/or modify it under the
 *  terms of the GNU General Public License as published by the Free Software
 *  Foundation, either version 3 of the License, or (at your option) any later
 *  version.
 *
 *  This software is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 *  more details.
 */

#ifndef DNA_FOUR_BIT_HPP
#define DNA_FOUR_BIT_HPP

#include <cstdint>
#include <iterator>
#include <vector>

enum base_in_byte { left, right };

extern char dna_four_bit_decoding[16];

template <typename uint_type> constexpr
uint_type
get_nibble(const uint_type x, const size_t offset) {
  return (x >> (4*offset)) & 15ul;
}

template <typename uint_type> constexpr
char
decode_dna_four_bit(const uint_type x,
                    const size_t offset) {
  return dna_four_bit_decoding[get_nibble(x, offset)];
}

template<class InputItr, class OutputIt>
OutputIt
decode_dna_four_bit(InputItr first, InputItr last, OutputIt d_first) {
  // ADS: assume destination has enough space
  while (first != last) {
    for (size_t offset = 0; offset < 16; ++offset)
      *d_first++ = decode_dna_four_bit(*first, offset);
    ++first;
  }
  // if original sequence length is odd and encoding not padded at the front,
  // then the final element in dest will be 'Z'
  return d_first;
}

template<class InCtr, class OutCtr>
void
decode_dna_four_bit(const InCtr &source, OutCtr &dest) {
  // expand out the bytes as pairs (do this backwards in case source == dest)
  const size_t source_size = source.size();
  dest.resize(16*source_size);
  size_t i = source_size;
  size_t j = dest.size();
  while (i > 0) {
    dest[--j] = source[--i];
    dest[--j] = source[i];
  }
  for (i = 0; i < dest.size(); i += 16) {
    for (size_t offset = 0; offset < 16; ++offset)
      dest[i + offset] = decode_dna_four_bit(dest[i], offset);
  }
}

extern uint8_t dna_four_bit_encoding[128];
template <typename uint_type> constexpr
size_t
encode_dna_four_bit(const uint_type x,
                    const size_t offset) {
  return (static_cast<size_t>(
           dna_four_bit_encoding[static_cast<unsigned>(x)])
         ) << (4*offset);
}

template<class InputItr, class OutputIt>
OutputIt
encode_dna_four_bit(InputItr first, InputItr last, OutputIt d_first) {
  while (first != last) {
    *d_first = 0;
    for (size_t i = 0; i < 16 && first != last; ++i)
      *d_first |= encode_dna_four_bit(std::move(*first++), i);
    ++d_first;
  }
  return d_first;
}

// GS: intended to be used as pointer to 4-bit encoding of DNA within a vector
// of size_t values
struct genome_four_bit_itr {
  genome_four_bit_itr(const std::vector<size_t>::const_iterator itr_,
                      const int off_ = 0) : itr(itr_), offset(off_) {}

  size_t operator*() const {
    return (*itr >> (offset << 2)) & 15ul;
  }
  genome_four_bit_itr& operator++() {
    offset = (offset + 1) & 15ul;
    itr += (offset == 0);
    return *this;
  }
  genome_four_bit_itr operator++(int) {
    genome_four_bit_itr tmp(*this);
    offset = (offset + 1) & 15ul;
    itr += (offset == 0);
    return tmp;
  }
  genome_four_bit_itr& operator--() {
    itr -= (offset == 0);

    offset = (offset - 1) & 15ul;
    return *this;
  }
  genome_four_bit_itr operator--(int) {
    genome_four_bit_itr tmp(*this);
    itr -= (offset == 0);
    offset = (offset - 1) & 15ul;
    return tmp;
  }
  genome_four_bit_itr operator+(const size_t step) const {
    // whether the sum of offsets is >= 16
    const bool shift_one_pos =
      (((offset + (static_cast<int>(step) & 15)) & 16) >> 4);

    const int new_offset = (offset + step) & 15;
    return genome_four_bit_itr(itr + step/16 + shift_one_pos,
                               new_offset);
  }
  bool operator!=(const genome_four_bit_itr &rhs) const {
    return itr != rhs.itr || offset != rhs.offset;
  }
  std::vector<size_t>::const_iterator itr;
  int offset;
};

#endif
