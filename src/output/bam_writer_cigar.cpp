// LLmap — BAM/SAM CIGAR string utilities.

#include "output/bam_writer.h"

#include <string>
#include <string_view>

namespace llmap::output::cigar {

std::string GenerateSimpleCigar(std::uint64_t query_len,
                                 std::uint64_t ref_start,
                                 std::uint64_t ref_end,
                                 std::uint32_t nm) {
    // Simple heuristic: if aligned length matches query and no edits, all M
    std::uint64_t ref_len = ref_end - ref_start;

    if (nm == 0 && ref_len == query_len) {
        return std::to_string(query_len) + "M";
    }

    // If edit distance exists, approximate with matches and mismatches
    if (ref_len == query_len) {
        // Same length: all differences are substitutions
        return std::to_string(query_len) + "M";
    }

    // Different lengths: needs indels
    if (ref_len > query_len) {
        // Deletions needed
        std::uint64_t del_len = ref_len - query_len;
        std::uint64_t match_before = query_len / 2;
        std::uint64_t match_after = query_len - match_before;
        return std::to_string(match_before) + "M" +
               std::to_string(del_len) + "D" +
               std::to_string(match_after) + "M";
    } else {
        // Insertions needed
        std::uint64_t ins_len = query_len - ref_len;
        std::uint64_t match_before = ref_len / 2;
        std::uint64_t match_after = ref_len - match_before;
        return std::to_string(match_before) + "M" +
               std::to_string(ins_len) + "I" +
               std::to_string(match_after) + "M";
    }
}

std::pair<std::uint64_t, std::uint64_t> CigarStats(std::string_view cigar) {
    std::uint64_t query_consumed = 0;
    std::uint64_t ref_consumed = 0;

    std::uint64_t num = 0;
    for (char c : cigar) {
        if (c >= '0' && c <= '9') {
            num = num * 10 + (c - '0');
        } else {
            switch (c) {
                case 'M':
                case '=':
                case 'X':
                    query_consumed += num;
                    ref_consumed += num;
                    break;
                case 'I':
                case 'S':
                    query_consumed += num;
                    break;
                case 'D':
                case 'N':
                    ref_consumed += num;
                    break;
                case 'H':
                case 'P':
                    break;
            }
            num = 0;
        }
    }

    return {query_consumed, ref_consumed};
}

bool ValidateCigar(std::string_view cigar) {
    if (cigar.empty() || cigar == "*") return true;

    bool expect_num = true;
    bool has_num = false;

    for (char c : cigar) {
        if (c >= '0' && c <= '9') {
            if (!expect_num && !has_num) return false;
            expect_num = true;
            has_num = true;
        } else if (c == 'M' || c == 'I' || c == 'D' || c == 'N' ||
                   c == 'S' || c == 'H' || c == 'P' || c == '=' || c == 'X') {
            if (!has_num) return false;
            has_num = false;
            expect_num = true;
        } else {
            return false;  // Invalid character
        }
    }

    return !has_num;  // Should end with an operation, not a number
}

}  // namespace llmap::output::cigar
