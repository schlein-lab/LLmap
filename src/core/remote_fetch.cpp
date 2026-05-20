// LLmap — remote reference/anchor fetching implementation.

#include "core/remote_fetch.h"

#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <string>
#include <vector>

namespace llmap::core {

namespace {

bool StartsWith(std::string_view s, std::string_view p) {
    return s.size() >= p.size() && s.compare(0, p.size(), p) == 0;
}

// Run argv via fork/exec (no shell -> no injection). Returns true on exit 0.
bool RunExec(const std::vector<std::string>& argv) {
    std::vector<char*> cargv;
    cargv.reserve(argv.size() + 1);
    for (const auto& a : argv) cargv.push_back(const_cast<char*>(a.c_str()));
    cargv.push_back(nullptr);

    pid_t pid = fork();
    if (pid < 0) return false;
    if (pid == 0) {
        // child: silence stdout/stderr
        ::freopen("/dev/null", "w", stdout);
        ::freopen("/dev/null", "w", stderr);
        ::execvp(cargv[0], cargv.data());
        _exit(127);  // exec failed (e.g. downloader not installed)
    }
    int status = 0;
    if (::waitpid(pid, &status, 0) < 0) return false;
    return WIFEXITED(status) && WEXITSTATUS(status) == 0;
}

// Download url -> dest using curl, then wget. Returns true on success.
bool Download(const std::string& url, const std::string& dest) {
    if (RunExec({"curl", "-fsSL", "-o", dest, url})) return true;
    if (RunExec({"wget", "-q", "-O", dest, url})) return true;
    return false;
}

bool FileNonEmpty(const std::string& path) {
    struct stat st {};
    return ::stat(path.c_str(), &st) == 0 && st.st_size > 0;
}

std::string CacheDir(std::string_view requested) {
    if (!requested.empty()) return std::string(requested);
    if (const char* e = std::getenv("LLMAP_REMOTE_CACHE")) return e;
    std::string base = "/tmp";
    if (const char* t = std::getenv("TMPDIR")) base = t;
    return base + "/llmap_remote_cache";
}

std::string Basename(std::string_view url) {
    auto slash = url.find_last_of('/');
    std::string b = (slash == std::string_view::npos)
                        ? std::string(url)
                        : std::string(url.substr(slash + 1));
    // strip query string
    auto q = b.find('?');
    if (q != std::string::npos) b = b.substr(0, q);
    if (b.empty()) b = "remote";
    return b;
}

}  // namespace

bool IsRemotePath(std::string_view path) {
    return StartsWith(path, "http://") || StartsWith(path, "https://") ||
           StartsWith(path, "s3://");
}

std::string S3ToHttps(std::string_view url) {
    constexpr std::string_view kPrefix = "s3://";
    if (!StartsWith(url, kPrefix)) return std::string(url);
    std::string_view rest = url.substr(kPrefix.size());
    auto slash = rest.find('/');
    if (slash == std::string_view::npos) {
        return "https://" + std::string(rest) + ".s3.amazonaws.com/";
    }
    std::string_view bucket = rest.substr(0, slash);
    std::string_view key = rest.substr(slash + 1);
    return "https://" + std::string(bucket) + ".s3.amazonaws.com/" +
           std::string(key);
}

std::optional<std::string> FetchToCache(std::string_view url,
                                        std::string_view cache_dir,
                                        bool also_index) {
    const std::string http = IsRemotePath(url) ? S3ToHttps(url) : std::string(url);
    const std::string dir = CacheDir(cache_dir);
    std::error_code ec;
    std::filesystem::create_directories(dir, ec);

    const std::uint64_t h = std::hash<std::string>{}(http);
    char hashbuf[17];
    std::snprintf(hashbuf, sizeof(hashbuf), "%016llx",
                  static_cast<unsigned long long>(h));
    const std::string local = dir + "/" + hashbuf + "_" + Basename(http);

    if (!FileNonEmpty(local)) {
        if (!Download(http, local)) return std::nullopt;
    }
    if (also_index) {
        // best-effort: ignore failures (region queries just won't work)
        const std::string fai = local + ".fai";
        const std::string gzi = local + ".gzi";
        if (!FileNonEmpty(fai)) Download(http + ".fai", fai);
        if (!FileNonEmpty(gzi)) Download(http + ".gzi", gzi);
    }
    return local;
}

}  // namespace llmap::core
