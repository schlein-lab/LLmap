#include "classical/wfa2_aligner.h"

#include <chrono>
#include <sstream>

namespace llmap::classical {

char CigarElement::OpChar() const {
    switch (op) {
        case CigarOp::Match:     return 'M';
        case CigarOp::Insertion: return 'I';
        case CigarOp::Deletion:  return 'D';
        case CigarOp::SoftClip:  return 'S';
        case CigarOp::Equal:     return '=';
        case CigarOp::Diff:      return 'X';
        default: return '?';
    }
}

std::string CigarElement::ToString() const {
    return std::to_string(length) + OpChar();
}

std::string WFA2Result::CigarString() const {
    std::string result;
    for (const auto& elem : cigar) {
        result += elem.ToString();
    }
    return result;
}

size_t WFA2Result::AlignedLength() const {
    size_t len = 0;
    for (const auto& elem : cigar) {
        if (elem.op != CigarOp::Insertion && elem.op != CigarOp::SoftClip) {
            len += elem.length;
        }
    }
    return len;
}

float ComputeIdentity(const std::vector<CigarElement>& cigar) {
    size_t matches = 0;
    size_t aligned = 0;

    for (const auto& elem : cigar) {
        switch (elem.op) {
            case CigarOp::Match:
            case CigarOp::Equal:
                matches += elem.length;
                aligned += elem.length;
                break;
            case CigarOp::Diff:
                aligned += elem.length;
                break;
            case CigarOp::Insertion:
            case CigarOp::Deletion:
                aligned += elem.length;
                break;
            default:
                break;
        }
    }

    return aligned > 0 ? static_cast<float>(matches) / static_cast<float>(aligned) : 0.0f;
}

std::optional<std::vector<CigarElement>> ParseCigar(std::string_view cigar_str) {
    std::vector<CigarElement> result;
    size_t i = 0;

    while (i < cigar_str.size()) {
        // Parse length
        size_t start = i;
        while (i < cigar_str.size() && std::isdigit(cigar_str[i])) {
            ++i;
        }

        if (i == start || i >= cigar_str.size()) {
            return std::nullopt;
        }

        uint32_t length = 0;
        for (size_t j = start; j < i; ++j) {
            length = length * 10 + (cigar_str[j] - '0');
        }

        // Parse operation
        CigarOp op;
        switch (cigar_str[i]) {
            case 'M': op = CigarOp::Match; break;
            case 'I': op = CigarOp::Insertion; break;
            case 'D': op = CigarOp::Deletion; break;
            case 'S': op = CigarOp::SoftClip; break;
            case '=': op = CigarOp::Equal; break;
            case 'X': op = CigarOp::Diff; break;
            default: return std::nullopt;
        }

        result.push_back({op, length});
        ++i;
    }

    return result;
}

bool IsWFA2Available() {
#ifdef LLMAP_HAS_WFA2
    return true;
#else
    return false;
#endif
}

}  // namespace llmap::classical
