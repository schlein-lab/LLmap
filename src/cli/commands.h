// LLmap — CLI command declarations.
//
// Each command is implemented in its own file (cmd_*.cpp).

#pragma once

namespace llmap::cli {

int run_allpair(int argc, char** argv);
int run_generate_synth(int argc, char** argv);
int run_validate_real(int argc, char** argv);
int run_index(int argc, char** argv);
int run_align(int argc, char** argv);
int run_sc_paralog_matrix(int argc, char** argv);
int run_sc_qc_report(int argc, char** argv);
int run_check(int argc, char** argv);
int run_annotate_ref(int argc, char** argv);

}  // namespace llmap::cli
