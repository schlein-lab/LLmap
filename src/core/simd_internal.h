#pragma once

#include <cstdint>

namespace llmap::core::simd {

// Lookup table for base encoding: ASCII -> 2-bit (A=0, C=1, G=2, T=3)
extern const uint8_t kBaseEncode[256];

// Validity table: 1 for valid ACGT bases, 0 for everything else
extern const uint8_t kBaseValid[256];

}  // namespace llmap::core::simd
