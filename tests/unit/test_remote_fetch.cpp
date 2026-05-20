// LLmap — unit tests for remote path detection and S3->HTTPS translation.
// Network-free: exercises only the pure helpers.

#include "core/remote_fetch.h"

#include <gtest/gtest.h>

using namespace llmap::core;

TEST(RemoteFetch, IsRemotePath) {
    EXPECT_TRUE(IsRemotePath("http://example.org/x.fa"));
    EXPECT_TRUE(IsRemotePath("https://example.org/x.fa"));
    EXPECT_TRUE(IsRemotePath("s3://bucket/key.fa.gz"));
    EXPECT_FALSE(IsRemotePath("/local/path/x.fa"));
    EXPECT_FALSE(IsRemotePath("relative/x.fa"));
    EXPECT_FALSE(IsRemotePath("ftp://example.org/x.fa"));  // not supported
    EXPECT_FALSE(IsRemotePath(""));
}

TEST(RemoteFetch, S3ToHttps) {
    EXPECT_EQ(S3ToHttps("s3://human-pangenomics/working/HPRC/HG02027/x.fa.gz"),
              "https://human-pangenomics.s3.amazonaws.com/working/HPRC/HG02027/x.fa.gz");
    EXPECT_EQ(S3ToHttps("s3://bucket/key"),
              "https://bucket.s3.amazonaws.com/key");
    // bucket only
    EXPECT_EQ(S3ToHttps("s3://bucket"),
              "https://bucket.s3.amazonaws.com/");
    // non-s3 passthrough
    EXPECT_EQ(S3ToHttps("https://example.org/x.fa"),
              "https://example.org/x.fa");
    EXPECT_EQ(S3ToHttps("/local/x.fa"), "/local/x.fa");
}
