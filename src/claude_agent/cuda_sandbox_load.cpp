#include "cuda_sandbox.h"

#include <atomic>
#include <chrono>
#include <fstream>
#include <map>
#include <sstream>

namespace llmap::claude_agent {

namespace {

std::string GenerateModuleId() {
    static std::atomic<uint64_t> counter{0};
    auto now = std::chrono::system_clock::now().time_since_epoch().count();
    char buf[32];
    snprintf(buf, sizeof(buf), "mod_%lu_%lu",
             static_cast<unsigned long>(now),
             static_cast<unsigned long>(++counter));
    return std::string(buf);
}

}  // namespace

struct CudaLoader::Impl {
    KernelBudget budget;
    std::map<std::string, std::vector<LoadedKernel>> modules;

    explicit Impl(KernelBudget b) : budget(std::move(b)) {}
};

CudaLoader::CudaLoader(KernelBudget budget)
    : impl_(std::make_unique<Impl>(std::move(budget))) {}

CudaLoader::~CudaLoader() = default;

CudaLoader::CudaLoader(CudaLoader&&) noexcept = default;
CudaLoader& CudaLoader::operator=(CudaLoader&&) noexcept = default;

LoadResult CudaLoader::Load(const std::filesystem::path& cubin_path) {
    auto start = std::chrono::steady_clock::now();
    LoadResult result;

    if (!std::filesystem::exists(cubin_path)) {
        result.status = LoadStatus::FileNotFound;
        return result;
    }

#ifdef LLMAP_HAS_CUDA
    CUmodule module;
    CUresult res = cuModuleLoad(&module, cubin_path.string().c_str());
    if (res != CUDA_SUCCESS) {
        result.status = LoadStatus::InvalidModule;
        return result;
    }

    std::vector<LoadedKernel> kernels;

    int func_count;
    res = cuModuleGetFunctionCount(&func_count, module);
    if (res == CUDA_SUCCESS) {
        for (int i = 0; i < func_count; ++i) {
            CUfunction func;
            res = cuModuleGetFunction(&func, module, kernel_name.c_str());
            if (res == CUDA_SUCCESS) {
                LoadedKernel k;
                k.function_ptr = reinterpret_cast<void*>(func);

                int regs, shared, max_threads;
                cuFuncGetAttribute(&regs, CU_FUNC_ATTRIBUTE_NUM_REGS, func);
                cuFuncGetAttribute(&shared, CU_FUNC_ATTRIBUTE_SHARED_SIZE_BYTES, func);
                cuFuncGetAttribute(&max_threads, CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK, func);

                k.registers_per_thread = regs;
                k.shared_memory_bytes = shared;
                k.max_threads_per_block = max_threads;

                if (k.registers_per_thread > impl_->budget.max_registers_per_thread ||
                    k.shared_memory_bytes > impl_->budget.max_shared_memory_bytes) {
                    result.status = LoadStatus::BudgetExceeded;
                    cuModuleUnload(module);
                    return result;
                }

                kernels.push_back(std::move(k));
            }
        }
    }

    std::string module_id = GenerateModuleId();
    impl_->modules[module_id] = std::move(kernels);

    result.status = LoadStatus::Success;
    result.module_id = module_id;
    result.kernels = impl_->modules[module_id];
#else
    auto module_id = GenerateModuleId();
    impl_->modules[module_id] = {};

    result.status = LoadStatus::Success;
    result.module_id = module_id;
#endif

    auto end = std::chrono::steady_clock::now();
    result.load_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    return result;
}

void CudaLoader::Unload(std::string_view module_id) {
    auto it = impl_->modules.find(std::string(module_id));
    if (it != impl_->modules.end()) {
        impl_->modules.erase(it);
    }
}

void CudaLoader::UnloadAll() {
    impl_->modules.clear();
}

std::optional<LoadedKernel> CudaLoader::GetKernel(std::string_view module_id,
                                                   std::string_view kernel_name) const {
    auto it = impl_->modules.find(std::string(module_id));
    if (it == impl_->modules.end()) return std::nullopt;

    for (const auto& kernel : it->second) {
        if (kernel.name == kernel_name) {
            return kernel;
        }
    }
    return std::nullopt;
}

const KernelBudget& CudaLoader::GetBudget() const {
    return impl_->budget;
}

size_t CudaLoader::LoadedModuleCount() const {
    return impl_->modules.size();
}

std::vector<std::string> CudaLoader::LoadedModuleIds() const {
    std::vector<std::string> ids;
    ids.reserve(impl_->modules.size());
    for (const auto& [id, _] : impl_->modules) {
        ids.push_back(id);
    }
    return ids;
}

struct CudaSandbox::Impl {
    Config config;
    CudaAnalyzer analyzer;
    CudaCompiler compiler;
    CudaLoader loader;
    std::vector<AuditEntry> audit_log;
    std::map<std::string, std::string> kernel_to_module;

    explicit Impl(Config cfg)
        : config(std::move(cfg))
        , analyzer(config.analyzer)
        , compiler(config.compiler)
        , loader(config.budget)
    {}
};

CudaSandbox::CudaSandbox(Config config)
    : impl_(std::make_unique<Impl>(std::move(config))) {}

CudaSandbox::~CudaSandbox() = default;

CudaSandbox::CudaSandbox(CudaSandbox&&) noexcept = default;
CudaSandbox& CudaSandbox::operator=(CudaSandbox&&) noexcept = default;

CudaSandbox::ExecutionResult CudaSandbox::CompileAndLoad(std::string_view cuda_source,
                                                          std::string_view kernel_name) {
    ExecutionResult result;
    std::string source_hash = ComputeSourceHash(cuda_source);

    AuditEntry entry{
        .timestamp = std::chrono::system_clock::now(),
        .session_id = impl_->config.session_id,
        .action = "compile_and_load",
        .allowed = false,
        .source_hash = source_hash,
        .details = "kernel=" + std::string(kernel_name)
    };

    result.analysis = impl_->analyzer.Analyze(cuda_source);

    if (impl_->config.require_analysis_pass && !result.analysis.safe) {
        result.success = false;
        result.error_message = "Static analysis failed: " +
            std::to_string(result.analysis.violations.size()) + " violations";
        entry.details += "; analysis_failed";
        WriteAuditEntry(entry);
        return result;
    }

    result.compile = impl_->compiler.Compile(cuda_source, std::string(kernel_name));

    if (result.compile.status != CompileStatus::Success) {
        result.success = false;
        result.error_message = "Compilation failed: " +
            std::string(CompileStatusName(result.compile.status));
        entry.details += "; compile_failed=" + result.compile.stderr_log;
        WriteAuditEntry(entry);
        return result;
    }

    result.load = impl_->loader.Load(*result.compile.output_path);

    if (result.load.status != LoadStatus::Success) {
        result.success = false;
        result.error_message = "Load failed: " +
            std::string(LoadStatusName(result.load.status));
        entry.details += "; load_failed";
        WriteAuditEntry(entry);
        return result;
    }

    impl_->kernel_to_module[std::string(kernel_name)] = result.load.module_id;

    result.success = true;
    entry.allowed = true;
    entry.details += "; success";
    WriteAuditEntry(entry);

    return result;
}

std::optional<LoadedKernel> CudaSandbox::GetKernel(std::string_view kernel_name) const {
    auto it = impl_->kernel_to_module.find(std::string(kernel_name));
    if (it == impl_->kernel_to_module.end()) return std::nullopt;
    return impl_->loader.GetKernel(it->second, kernel_name);
}

void CudaSandbox::WriteAuditEntry(const AuditEntry& entry) {
    impl_->audit_log.push_back(entry);

    if (!impl_->config.audit_log_path.empty()) {
        std::ofstream out(impl_->config.audit_log_path, std::ios::app);
        if (out) {
            auto time_t = std::chrono::system_clock::to_time_t(entry.timestamp);
            out << std::ctime(&time_t);
            out << "session=" << entry.session_id << "\n";
            out << "action=" << entry.action << "\n";
            out << "allowed=" << (entry.allowed ? "yes" : "no") << "\n";
            out << "hash=" << entry.source_hash << "\n";
            out << "details=" << entry.details << "\n";
            out << "---\n";
        }
    }
}

std::vector<AuditEntry> CudaSandbox::ReadAuditLog() const {
    return impl_->audit_log;
}

void CudaSandbox::ClearAuditLog() {
    impl_->audit_log.clear();
}

const CudaSandbox::Config& CudaSandbox::GetConfig() const {
    return impl_->config;
}

}  // namespace llmap::claude_agent
