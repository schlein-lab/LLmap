// LLmap — Stage2Pipeline internal implementation details.
// Split across: stage2_pipeline_core.cpp, stage2_pipeline_init.cpp, stage2_pipeline_em.cpp

#pragma once

#include "reference_collapse/stage2_pipeline.h"
#include "reference_collapse/collapse_check.h"
#include "reference_collapse/em_iterator.h"
#include "reference_collapse/member_propagation.h"
#include "reference_collapse/reference_index.h"
#include "reference_collapse/refinement.h"
#include "core/wave_state.h"

#include <chrono>

namespace llmap {

// Internal state for Stage 2 pipeline
struct Stage2Pipeline::InternalState {
    std::unique_ptr<ReferenceIndex> ref_index;
    std::unique_ptr<WaveState> wave_state;
    std::unique_ptr<Stage2Result> result;

    EmIterator em_iterator;
    CollapseChecker collapse_checker;
    Refinement refinement;
    MemberPropagation member_propagation;

    ChildIndex l0_to_l1_index;
    ChildIndex l1_to_l2_index;

    Stage2Stats stats;
    std::string last_error;

    std::chrono::steady_clock::time_point start_time;
};

}  // namespace llmap
