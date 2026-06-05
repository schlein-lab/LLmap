// LLmap — Internal types for memory-mapped FASTA reader.
// This file is internal to the mmap_fasta implementation.

#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace llmap::io {

// Internal sequence entry with full metadata
struct SequenceEntry {
    std::string name;           // Owned copy of name
    size_t name_offset;         // Offset of '>' in file
    size_t data_offset;         // Offset where sequence data starts
    size_t data_end;            // Offset where sequence data ends
    size_t length;              // Actual nucleotide count (excluding newlines)
};

// PIMPL implementation class
class MmapFastaImpl {
public:
    int fd = -1;
    void* mapped = nullptr;
    size_t file_size = 0;
    std::string last_error;

    // When the input file was gz-compressed we decompress the entire
    // payload into this heap buffer at open-time and point `mapped` at
    // its data. The destructor must delete[] the buffer instead of
    // calling munmap, and madvise calls become no-ops.
    bool   owned_heap = false;
    char*  heap_buffer = nullptr;

    std::vector<SequenceEntry> sequences;
    std::unordered_map<std::string, size_t> name_to_index;

    ~MmapFastaImpl();

    const char* Data() const;
    bool BuildIndex();
    std::string ExtractSequence(size_t data_offset, size_t data_end) const;
    std::string ExtractSubsequence(size_t data_offset, size_t data_end,
                                   size_t start, size_t len) const;
};

// Get system page size
size_t GetPageSize();

}  // namespace llmap::io
