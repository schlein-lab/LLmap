#pragma once

#include "agent_types.h"

#include <filesystem>
#include <string>

namespace llmap::claude_agent {

std::string SerializeBiologyPrior(const BiologyPrior& prior);
BiologyPrior DeserializeBiologyPrior(std::string_view json);

bool WriteBiologyPrior(
    const BiologyPrior& prior,
    const std::filesystem::path& output_path
);

std::optional<BiologyPrior> ReadBiologyPrior(
    const std::filesystem::path& input_path
);

std::string SerializeSampleParams(const SampleParams& params);
SampleParams DeserializeSampleParams(std::string_view json);

bool WriteSampleParams(
    const SampleParams& params,
    const std::filesystem::path& output_path
);

std::optional<SampleParams> ReadSampleParams(
    const std::filesystem::path& input_path
);

}  // namespace llmap::claude_agent
