#include "annot/annot_types.h"

namespace llmap::annot {

void ParamOverride::OverlayOver(const ParamOverride& base) {
    if (!k                       && base.k)                       k = base.k;
    if (!w                       && base.w)                       w = base.w;
    if (!max_occ                 && base.max_occ)                 max_occ = base.max_occ;
    if (!lambda_scale            && base.lambda_scale)            lambda_scale = base.lambda_scale;
    if (!identity_threshold      && base.identity_threshold)      identity_threshold = base.identity_threshold;
    if (!anchor_weight_scale     && base.anchor_weight_scale)     anchor_weight_scale = base.anchor_weight_scale;
    if (!report_multi_position   && base.report_multi_position)   report_multi_position = base.report_multi_position;
    if (!require_psv_disambig    && base.require_psv_disambig)    require_psv_disambig = base.require_psv_disambig;
    if (!allow_high_mismatch     && base.allow_high_mismatch)     allow_high_mismatch = base.allow_high_mismatch;
    if (!require_llm_at_runtime  && base.require_llm_at_runtime)  require_llm_at_runtime = base.require_llm_at_runtime;
}

}  // namespace llmap::annot
