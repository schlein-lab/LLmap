#include "cuda_sandbox.h"

#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>

#ifdef __linux__
#include <sys/wait.h>
#include <unistd.h>
#endif

namespace llmap::claude_agent {

namespace {

std::string ExecuteCommand(const std::string& cmd, int& exit_code, size_t timeout_seconds) {
    std::string result;
    exit_code = -1;

#ifdef __linux__
    std::array<char, 4096> buffer;
    std::string full_cmd = "timeout " + std::to_string(timeout_seconds) + " " + cmd + " 2>&1";
    FILE* pipe = popen(full_cmd.c_str(), "r");
    if (!pipe) return "Failed to execute command";

    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) {
        result += buffer.data();
    }

    int status = pclose(pipe);
    if (WIFEXITED(status)) {
        exit_code = WEXITSTATUS(status);
    }
#else
    result = "Command execution not supported on this platform";
#endif

    return result;
}

bool FileExists(const std::filesystem::path& path) {
    return std::filesystem::exists(path);
}

void WriteSourceFile(const std::filesystem::path& path, std::string_view content) {
    std::ofstream out(path);
    out << content;
}

std::string BuildNvccCommand(const CompileConfig& config,
                             const std::filesystem::path& source,
                             const std::filesystem::path& output) {
    std::ostringstream cmd;

    cmd << config.nvcc_path.string();
    cmd << " -arch=" << config.arch;
    cmd << " -cubin";
    cmd << " -O2";
    cmd << " --disable-warnings";
    cmd << " -Xcompiler \"-fPIC\"";

    for (const auto& flag : config.extra_flags) {
        cmd << " " << flag;
    }

    cmd << " -o " << output.string();
    cmd << " " << source.string();

    return cmd.str();
}

std::string BuildBubblewrapCommandImpl(const CompileConfig& config,
                                       const std::filesystem::path& source,
                                       const std::filesystem::path& output) {
    std::ostringstream cmd;

    cmd << "bwrap";

    cmd << " --ro-bind /usr /usr";
    cmd << " --ro-bind /lib /lib";
    cmd << " --ro-bind /lib64 /lib64";
    cmd << " --ro-bind /bin /bin";
    cmd << " --ro-bind /sbin /sbin";

    if (FileExists("/etc/alternatives")) {
        cmd << " --ro-bind /etc/alternatives /etc/alternatives";
    }

    if (!config.nvcc_path.empty() && FileExists(config.nvcc_path.parent_path())) {
        auto cuda_root = config.nvcc_path.parent_path().parent_path();
        cmd << " --ro-bind " << cuda_root.string() << " " << cuda_root.string();
    }

    cmd << " --bind " << config.sandbox_dir.string() << " " << config.sandbox_dir.string();
    cmd << " --bind " << config.output_dir.string() << " " << config.output_dir.string();

    cmd << " --dev /dev";
    cmd << " --proc /proc";
    cmd << " --tmpfs /tmp";

    cmd << " --unshare-net";
    cmd << " --unshare-pid";
    cmd << " --unshare-uts";
    cmd << " --unshare-ipc";

    cmd << " --die-with-parent";
    cmd << " --new-session";

    cmd << " --setenv PATH /usr/local/cuda/bin:/usr/bin:/bin";

    cmd << " -- ";
    cmd << BuildNvccCommand(config, source, output);

    return cmd.str();
}

}  // namespace

struct CudaCompiler::Impl {
    CompileConfig config;
    bool bwrap_available{false};

    explicit Impl(CompileConfig cfg) : config(std::move(cfg)) {
#ifdef __linux__
        int exit_code;
        ExecuteCommand("which bwrap", exit_code, 5);
        bwrap_available = (exit_code == 0);
#endif
        std::filesystem::create_directories(config.sandbox_dir);
        std::filesystem::create_directories(config.output_dir);
    }
};

CudaCompiler::CudaCompiler(CompileConfig config)
    : impl_(std::make_unique<Impl>(std::move(config))) {}

CudaCompiler::~CudaCompiler() = default;

CudaCompiler::CudaCompiler(CudaCompiler&&) noexcept = default;
CudaCompiler& CudaCompiler::operator=(CudaCompiler&&) noexcept = default;

CompileResult CudaCompiler::Compile(std::string_view cuda_source, std::string_view output_name) {
    auto start = std::chrono::steady_clock::now();

    auto source_path = impl_->config.sandbox_dir / (std::string(output_name) + ".cu");
    WriteSourceFile(source_path, cuda_source);

    auto result = CompileFile(source_path, output_name);

    std::filesystem::remove(source_path);

    return result;
}

CompileResult CudaCompiler::CompileFile(const std::filesystem::path& source_path,
                                        std::string_view output_name) {
    auto start = std::chrono::steady_clock::now();
    CompileResult result;

    if (!FileExists(source_path)) {
        result.status = CompileStatus::CompileError;
        result.stderr_log = "Source file does not exist: " + source_path.string();
        return result;
    }

    auto output_path = impl_->config.output_dir / (std::string(output_name) + ".cubin");

    std::string cmd;
    if (impl_->config.use_bubblewrap && impl_->bwrap_available) {
        cmd = BuildBubblewrapCommandImpl(impl_->config, source_path, output_path);
    } else {
        cmd = BuildNvccCommand(impl_->config, source_path, output_path);
    }

    int exit_code;
    std::string output = ExecuteCommand(cmd, exit_code, impl_->config.timeout_seconds);

    auto end = std::chrono::steady_clock::now();
    result.compile_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    if (exit_code == 124) {
        result.status = CompileStatus::Timeout;
        result.stderr_log = "Compilation timed out after " +
                           std::to_string(impl_->config.timeout_seconds) + " seconds";
        return result;
    }

    if (exit_code != 0) {
        result.status = CompileStatus::CompileError;
        result.stderr_log = output;
        return result;
    }

    if (!FileExists(output_path)) {
        result.status = CompileStatus::CompileError;
        result.stderr_log = "Compilation completed but output file not found";
        return result;
    }

    result.status = CompileStatus::Success;
    result.output_path = output_path;
    result.stdout_log = output;

    return result;
}

const CompileConfig& CudaCompiler::GetConfig() const {
    return impl_->config;
}

bool CudaCompiler::IsBubblewrapAvailable() const {
    return impl_->bwrap_available;
}

std::string CudaCompiler::BuildBubblewrapCommand(const std::filesystem::path& source,
                                                  const std::filesystem::path& output) const {
    return BuildBubblewrapCommandImpl(impl_->config, source, output);
}

}  // namespace llmap::claude_agent
