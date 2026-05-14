// LLmap — Internal utilities for cb_preservation module.
// Not part of the public API.

#pragma once

#include <string_view>
#include <tuple>

namespace llmap::singlecell::internal {

// Parse a single SAM tag field: "XX:T:VALUE"
// Returns (tag_name, type_char, value). Empty tag_name on parse failure.
std::tuple<std::string, char, std::string> ParseSamTag(std::string_view field);

// Check if character is valid DNA base (ACGT, case-insensitive)
bool IsValidDnaBase(char c) noexcept;

// Check if character is valid DNA base or N
bool IsValidDnaBaseOrN(char c) noexcept;

}  // namespace llmap::singlecell::internal
