// LLmap — junction_hunter: multi-k consensus + junction caller.
//
// For one read against one NAHR-pair, run the per-position consensus
// pass and emit the JunctionRecord with the final call.
//
// At each read position we have up to 5 (k, hash) emissions from the
// tiler. Each hash is looked up in the pair-index's per-k tables to
// get a LocusClass. The position's consensus is the LocusClass that is
// reported by at least `consensus_min` (default 3) of the 5 k-values;
// if no class meets that threshold, the position is `Ambiguous`.
//
// The read-level classification then runs the outside-in convergence
// check on the consensus sequence:
//   - LEFT half: consensus positions should be predominantly LcrUp,
//     and the LcrUp genomic positions should increase monotonically
//     in read order.
//   - RIGHT half: same for LcrDown.
//   - Interior consensus positions anywhere on the read → CanonicalInterior.
//   - All same-class consensus → CanonicalUp / CanonicalDown.
//   - Mixed without monotonicity → ChimeraArtifact.
//   - Excessive Ambiguous → ParalogAmbiguous.

#pragma once

#include "junction_hunter/junction_hunter_types.h"
#include "junction_hunter/pair_kmer_index.h"
#include "junction_hunter/read_tiler.h"

namespace llmap::junction_hunter {

/// Single-(read, pair) caller. Returns a fully populated JunctionRecord.
JunctionRecord CallJunction(std::string_view read_id,
                            const ReadTiling& tiling,
                            const PairKmerIndex& pair_index,
                            const NahrPair& pair,
                            const MultiKConfig& cfg = {});

}  // namespace llmap::junction_hunter
